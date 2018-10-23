#include <limits>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include "cmdline.h"
#include "vcf.h"
#include "pheno.h"
#include "statsutil.h"
#include "lsfit.h"


using std::size_t;


namespace {

struct Parameter
{
    std::string vcf;
    std::string pheno;
    std::string covar;
    std::string out;
} par ;

template<typename T>
std::vector<T> intersect(std::vector<T> a, std::vector<T> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    std::vector<T> c;

    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(c));

    return c;
}

template<typename T1, typename T2>
std::vector<T1> subset(const std::vector<T1> &vec, const std::vector<T2> &idx)
{
    std::vector<T1> out;

    out.reserve(idx.size());
    for (auto i : idx)
        out.push_back(vec[i]);

    return out;
}

template<typename T>
std::vector<T> subset(const std::vector<T> &vec, const std::vector<char> &mask)
{
    std::vector<T> out;

    out.reserve( std::count(mask.begin(), mask.end(), 1) );

    auto n = vec.size();
    for (size_t i = 0; i < n; ++i) {
        if ( mask[i] )
            out.push_back(vec[i]);
    }

    return out;
}

void filter_ind(const std::vector<int> &idx, Phenotype &pt)
{
    subset(pt.ind,idx).swap(pt.ind);

    if ( ! pt.env.empty() )
        subset(pt.env,idx).swap(pt.env);

    if ( ! pt.blk.empty() )
        subset(pt.blk,idx).swap(pt.blk);

    for (auto &e : pt.dat)
        subset(e,idx).swap(e);
}

void filter_ind(const std::vector<int> &idx, Covariate &ct)
{
    subset(ct.ind,idx).swap(ct.ind);

    for (auto &e : ct.dat)
        subset(e,idx).swap(e);
}

void merge(Genotype &gt, Phenotype &pt, Covariate &ct, std::vector<int> &gi)
{
    bool docovar = ! ct.phe.empty() && ! ct.ind.empty();

    bool domerge = gt.ind != pt.ind;
    if ( ! domerge && docovar )
        domerge = gt.ind != ct.ind;

    if ( ! domerge ) {
        gi.resize(gt.ind.size());
        std::iota(gi.begin(), gi.end(), 0);
        return;
    }

    std::cerr << "INFO: performing data intersection by individual...\n";

    auto ind = intersect(gt.ind, pt.ind);
    if ( docovar )
        ind = intersect(ind, ct.ind);

    std::vector<int> pi, ci;
    gi.clear();

    for (auto itr = pt.ind.begin(); itr != pt.ind.end(); ++itr) {
        if ( std::binary_search(ind.begin(), ind.end(), *itr) ) {
            pi.push_back( itr - pt.ind.begin() );
            if ( docovar ) {
                auto itr2 = std::find(ct.ind.begin(), ct.ind.end(), *itr);
                ci.push_back(itr2 - ct.ind.begin());
            }
            auto itr3 = std::find(gt.ind.begin(), gt.ind.end(), *itr);
            gi.push_back(itr3 - gt.ind.begin());
        }
    }

    filter_ind(pi, pt);

    if ( docovar )
        filter_ind(ci, ct);

    std::cerr << "INFO: there are " << ind.size() << " individuals after intersection\n";
}

void parse_envblk(const Phenotype &pt, std::vector< std::vector<double> > &ac, std::vector< std::vector<double> > &ic)
{
    std::vector< std::vector<double> > xenv, xblk;

    if ( ! pt.env.empty() ) {
        idummy1(factor(pt.env), xenv);
        ac.insert(ac.end(), xenv.begin(), xenv.end());
        ic.insert(ic.end(), xenv.begin(), xenv.end());
    }

    if ( ! pt.blk.empty() ) {
        idummy1(factor(pt.blk), xblk);
        if ( xenv.empty() )
            ac.insert(ac.end(), xblk.begin(), xblk.end());
        else {
            idummy3(factor(pt.env), xenv);
            for (auto &e : xenv) {
                for (auto v : xblk) {
                    std::transform(e.begin(), e.end(), v.begin(), v.begin(), std::multiplies<double>());
                    ac.push_back(v);
                }
            }
        }
    }
}

int assoc_glm(const Genotype &gt, const std::vector<int> &gi, const std::vector<double> &y,
              const std::vector< std::vector<double> > &ac,
              const std::vector< std::vector<double> > &ic,
              std::vector< std::vector<double> > &res)
{
    auto n = y.size();
    auto m = gt.dat.size();

    // modelP, modelR2, mainP, mainR2, intP, intR2
    res.assign(ic.empty() ? 2 : 6, std::vector<double>(m, std::numeric_limits<double>::quiet_NaN()));

    double sst = calc_css(y);

    std::vector<double> x0(n, 1);
    for (auto &v : ac)
        x0.insert(x0.end(), v.begin(), v.end());

    std::vector<double> b;
    double dfe0 = 0, sse0 = 0;
    lsfit(y, x0, b, dfe0, sse0);

    for (size_t j = 0; j < m; ++j) {
        std::vector<size_t> idx;
        std::vector< std::vector<double> > x1;

        if (gt.ploidy == 1) {
            std::vector<allele_t> g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = static_cast<size_t>(gi[i]);
                auto a = gt.dat[j][ii];
                if ( a ) {
                    g.push_back(a);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }
        else {
            std::vector< std::pair<allele_t,allele_t> > g;
            for (size_t i = 0; i < n; ++i) {
                auto ii = static_cast<size_t>(gi[i]);
                auto a = gt.dat[j][ii*2];
                auto b = gt.dat[j][ii*2+1];
                if ( a && b ) {
                    if (a > b)
                        std::swap(a, b);
                    g.emplace_back(a, b);
                    idx.push_back(i);
                }
            }
            idummy2(factor(g), x1);
        }

        if ( x1.empty() )
            continue;

        auto n2 = idx.size();

        auto jsst = sst;
        auto jdfe0 = dfe0;
        auto jsse0 = sse0;

        std::vector<double> x, y2;

        if (n2 == n) {
            y2 = y;
            x = x0;
        }
        else {
            y2 = subset(y, idx);

            jsst = calc_css(y2);

            x.assign(n2, 1);
            for (auto &e : ac) {
                auto z = subset(e, idx);
                x.insert(x.end(), z.begin(), z.end());
            }

            lsfit(y2, x, b, jdfe0, jsse0);
        }

        for (auto &e : x1)
            x.insert(x.end(), e.begin(), e.end());

        double jdfe1 = 0, jsse1 = 0;
        lsfit(y2, x, b, jdfe1, jsse1);

        auto dfx = jdfe0 - jdfe1;
        auto ssx = jsse0 - jsse1;

        if ( ic.empty() ) {
            if (dfx > 0 && ssx > 0 && jdfe1 > 0 && jsse1 > 0) {
                auto f = (ssx / dfx) / (jsse1 / jdfe1);
                res[0][j] = fpval(f, dfx, jdfe1);
                res[1][j] = ssx / jsst;
            }
        }
        else {
            for (auto &e : x1) {
                for (auto z : ic) {
                    if (n2 != n)
                        z = subset(z, idx);
                    std::transform(e.begin(), e.end(), z.begin(), z.begin(), std::multiplies<double>());
                    x.insert(x.end(), z.begin(), z.end());
                }
            }

            double jdfe2 = 0, jsse2 = 0;
            lsfit(y2, x, b, jdfe2, jsse2);

            auto dfm = jdfe0 - jdfe2;
            auto ssm = jsse0 - jsse2;

            if (dfm > 0 && ssm > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssm / dfm) / (jsse2 / jdfe2);
                res[0][j] = fpval(f, dfm, jdfe2);
                res[1][j] = ssm / jsst;
            }

            if (dfx > 0 && ssx > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssx / dfx) / (jsse2 / jdfe2);
                res[2][j] = fpval(f, dfx, jdfe2);
                res[3][j] = ssx / jsst;
            }

            auto dfxi = jdfe1 - jdfe2;
            auto ssxi = jsse1 - jsse2;

            if (dfxi > 0 && ssxi > 0 && jdfe2 > 0 && jsse2 > 0) {
                auto f = (ssxi / dfxi) / (jsse2 / jdfe2);
                res[4][j] = fpval(f, dfxi, jdfe2);
                res[5][j] = ssxi / jsst;
            }
        }
    }

    return 0;
}

} // namesapce


int glm_gwas(int argc, char *argv[])
{
    std::cerr << "GLM-GWAS 1.0 (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--pheno", "phenotype file", "");
    cmd.add("--covar", "covariate file", "");
    cmd.add("--out", "output file", "glm-gwas.out");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.pheno = cmd.get("--pheno");
    par.covar = cmd.get("--covar");
    par.out = cmd.get("--out");

    Genotype gt;
    Phenotype pt;
    Covariate ct;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::cerr << "INFO: reading phenotype file...\n";
    if (read_pheno(par.pheno, pt) != 0)
        return 1;
    std::cerr << "INFO: " << pt.ind.size() << " observations, " << pt.phe.size() << " traits\n";

    if ( ! par.covar.empty() ) {
        std::cerr << "INFO: reading covariate file...\n";
        if (read_covar(par.covar, ct) != 0)
            return 1;
        std::cerr << "INFO: " << ct.ind.size() << " individuals, " << ct.phe.size() << " covariates\n";
    }

    std::vector<int> gi;
    merge(gt, pt, ct, gi);

    if ( gi.empty() ) {
        std::cerr << "ERROR: no valid observations are found\n";
        return 1;
    }

    std::vector< std::vector<double> > ac, ic;
    parse_envblk(pt, ac, ic);

    ac.insert(ac.end(), ct.dat.begin(), ct.dat.end());

    auto nt = pt.phe.size();
    std::vector< std::vector<double> > res;

    for (size_t t = 0; t < nt; ++t) {
        std::vector<char> keep;
        for (auto e : pt.dat[t])
            keep.push_back( std::isfinite(e) );

        auto y = pt.dat[t];
        auto ac2 = ac;
        auto ic2 = ic;
        auto gi2 = gi;

        if (std::find(keep.begin(), keep.end(), 0) != keep.end()) {
            y = subset(y, keep);
            for (auto &e : ac2)
                e = subset(e, keep);
            for (auto &e : ic2)
                e = subset(e, keep);
            gi2 = subset(gi, keep);
        }

        std::vector< std::vector<double> > z;
        assoc_glm(gt, gi2, y, ac2, ic2, z);
        res.insert(res.end(), z.begin(), z.end());
    }

    std::ofstream ofs(par.out + ".pval");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pval" << "\n";
        return 1;
    }

    ofs << "Locus\tChromosome\tPosition";
    for (auto &e : pt.phe) {
        ofs << "\t" << e << ".modelP\t" << e << ".modelR2";
        if ( ! ic.empty() )
            ofs << "\t" << e << ".mainP\t" << e << ".mainR2\t" << e << ".intP\t" << e << ".intR2";
    }
    ofs << "\n";

    auto m = gt.dat.size();

    for (size_t j = 0; j < m; ++j) {
        ofs << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j];
        for (auto &v : res)
            ofs << "\t" << v[j];
        ofs << "\n";
    }

    return  0;
}
