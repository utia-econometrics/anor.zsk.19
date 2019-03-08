#ifndef EMIMULTI_H
#define EMIMULTI_H

#include <ostream>
#include "scenarios.h"
#include "lpsolvers.h"
#include "mmpcvar.h"
#include "distributions.h"

#define REDUCETRIVIALCONSTRAINTS

enum ewhattodo { ebasic, ewsscompute, efixeddemand, ednonstac, econtamination, emstac };

ewhattodo constexpr whattodo = emstac;




enum products {plates, cuts, profiles, brams,
              iron, numproducts, numsold=numproducts-1};

static const unsigned int maxfutures = 4;

static const unsigned int maxT = 4;

static const double varrho = 0.04;/// According to Smid/Zapletal 2017
static const double sigma = 0.15;/// According to Smid/Zapletal 2017
const double theta = 1e5; // penalization constant

extern double alphanstac;
extern double betanstac;


extern double alpha;
extern double beta;



// declarations of paste from problem params.xls

// (I-A)^{-1}
extern const double R[numproducts][numproducts];

// margins EUR/1000t
extern const double	m	[numproducts];

// production limits
extern const double	w	[numproducts];

// mean demand (1000t)
extern const double	Ed	[numproducts];

// log conveniences
extern const double	Ey[maxfutures];


// emisions (t/1000t final)
extern const double	hR	[numproducts];

// free alocated amounts (2017-2020)
extern const double	r	[maxT];

// spots
extern double p0;

// random vector	(d_cuts,d_profiles,d_brams,r)

extern std::vector<double> Exi1to3;

// root of variance matrix
extern std::vector<std::vector<double>> sqV;

extern const double fixedcost;

extern const double deadlyloss;


// atom of 3 point moment matching distribution
extern const double n3atom;

using emixi=std::vector<double>;

// declarations pasted from xls


// d[1], d[2], d[3], rho
extern const emixi xi0;

// end of declarations of paste from xls

inline double Exi3(const std::vector<double>& rdist, double stdev)
{
    assert(rdist.size());
    double s = 0;
    for(unsigned int i=0; i< rdist.size(); i++)
        s += exp(stdev*rdist[i]);
    return -log(s/rdist.size());
}


struct emiatom
{
    double d[numsold];
    double m[numsold];
    double p;
    double q[maxfutures];

    emiatom()
    {
        p=0;
        for(unsigned int i=0; i<numsold; i++)
            d[i] = m[i] = 0;
        for(unsigned int i=0; i<maxfutures; i++)
            q[i] = 0;
    }

    emiatom(emixi ex, double prevp,std::vector<double>& prevd)
    {
        assert(ex.size()==5);
        double mfactor = ex[4];
        for(unsigned int i=0; i<numsold; i++)
            m[i]=::m[i]*mfactor;
        double dr[numsold];
        for(unsigned int i=1; i<numsold; i++)
            dr[i]=ex[i-1];
        dr[0] = ex[0] * Ed[0] / Ed[1];
        bool minus = false;
        for(unsigned int i=0; i<numsold; i++)
           d[i]=alpha * prevd[i]
               + (1-alpha)*(Ed[i]+(dr[i]-Ed[i])*beta);
        p = prevp*exp(ex[3]);
        for(unsigned int i=0; i<maxfutures; i++)
            q[i]=p*exp(::Ey[i]);
    }
};


class emiproblem: public linearproblem<emiatom>
{
    friend class emitreesolutionstats;
    enum negoffsets {zno = 1, yno, cno, eno, sno };

    static unsigned int numfutures(unsigned int stage, unsigned int aT)
    {
        return std::min(maxfutures, aT - stage);
    }

    unsigned int numfutures(unsigned int stage) const
    {
        return numfutures(stage,T());
    }

    static std::vector<unsigned int> dims(unsigned int aT)
    {
        std::vector<unsigned int> r(aT+1);

        unsigned int i=0;
        for(; i<aT; i++)
            r[i] = numsold+numfutures(i,aT)+5; // s,e,c,y,z
        r[i] = 5;
        return r;
    }

public:
    emiproblem(
       bool z, // true/false - payments organized by time/measurability (apokryfy.lys)
       unsigned int aT, // horizon
       bool losspenalized, // deadly loss penalized by theta (gives probabily infeasible problems)
       bool zerofutures = false // true if zero number of futures included
         )
        : linearproblem<emiatom>(dims(aT)), fz(z), flp(losspenalized), fzf(zerofutures)
    {
        assert(aT<=maxT);
#ifdef REDUCETRIVIALCONSTRAINTS
        for(unsigned int i=0; i<numsold; i++)
        {
            trivialproductionconstraint[i] = true;
            if(R[i][i] != 1 ) // always should be 1, however, better safe than sorry
            {
                trivialproductionconstraint[i] = false;
                continue;
            }
            for(unsigned int j=0; j<numproducts; j++)
            {
                if(j!=i)
                    if(R[i][j] != 0)
                    {
                        trivialproductionconstraint[i] = false;
                        break;
                    }
            }
        }
        trivialproductionconstraint[numproducts-1] = false;
#else
        for(unsigned int j=0; j<numproducts; j++)
             trivialproductionconstraint[j] = false;
#endif
    }
    virtual void f(unsigned int stage,
                   const scenario<emiatom>& exi,
                     linearfunction_ptr& r) const
    {
        r.reset(new linearfunction(stagedim(stage)));

        assert(r->coefs.size());
        if(fz)
        {
            r->coefs[r->coefs.size()-zno ] = 1;
            if(stage==this->T())
                r->coefs[r->coefs.size()-cno ] = 1;
        }
        else
            r->coefs[r->coefs.size()- yno ] = 1;
    }

    virtual std::string varname(unsigned int stage, unsigned int i) const
    {
        bool laststage = stage==this->T();
        unsigned int np;
        if(laststage)
           np=0;
        else
        {
            switch(i)
            {
            case 0:
                return "plates";
            case 1:
                return "cuts";
            case 2:
                return "profiles";
            case 3:
                return "brams";
            }
            np = 4;
        }
        if(i < np+numfutures(stage))
        {
            static char n[3];
            n[0] = 'f';
            n[1] = '1'+(i-np);
            n[2] = 0;
            return n;
        }
        switch(this->stagedim(stage) -i)
        {
        case sno:
            return "s";
        case eno:
            return "e";
        case cno:
            return "c";
        case yno:
            return "y";
        case zno:
            return "z";
        }
        assert(1);
    }

    virtual void constraints(
                unsigned int stage,
                const scenario<emiatom>& xi,
                varinfo_list_ptr& vars,
                constraint_list_ptr<linearconstraint>& constraints
                ) const
    {
        bool laststage = stage==this->T();

        vars.reset(new varinfo_list);
        if(!laststage)
        {
            vars->push_back(varinfo(0,xi[stage].d[plates]));
            vars->push_back(varinfo(0,xi[stage].d[cuts]));
            vars->push_back(varinfo(0,xi[stage].d[profiles]));
            vars->push_back(varinfo(0,xi[stage].d[brams]));
            for(unsigned int i=0; i<numsold; i++)
            {
                if(trivialproductionconstraint[i]
                     && (*vars)[i].h > w[i])
                    (*vars)[i].h = w[i];
            }
            for(unsigned int i=0; i<numfutures(stage); i++)
                vars->push_back(fzf ? varinfo(0,0) : varinfo(varinfo::Rplus));
        }

        double rb = stage ? r[stage-1] :0;
        vars->push_back(varinfo( -rb,inf)); // s //laststage ? 0 :
        if(laststage)
            vars->push_back(varinfo(0,0)); //e
        else
            vars->push_back(varinfo(varinfo::Rplus)); //e
        vars->push_back(varinfo(varinfo::Rplus)); // c
        vars->push_back(varinfo(varinfo::R)); // z
        vars->push_back(varinfo(varinfo::R)); // y

        constraints.reset(new linearconstraint_list);

        if(!laststage)
        {
            for(unsigned int i=0; i<numproducts; i++)
            {
                if(!trivialproductionconstraint[i])
                {
                    constraints->push_back(linearconstraint(
                                                dimupto(stage),
                                                linearconstraint::leq,
                                                 w[i]));
                    for(unsigned int j=0; j<numsold; j++)
                    {
                        assert(stageoffset(stage)+j < constraints->rbegin()->lhs.size());
                        constraints->rbegin()->lhs[stageoffset(stage)+j]=R[i][j];
                    }
                }
            }
        }

        bool dlpenalized = flp && laststage;

        constraints->push_back(linearconstraint(
                                    dimupto(stage),
                                    linearconstraint::eq,
                                     rb));// e

        constraints->push_back(linearconstraint(
                                    dimupto(stage),
                                    linearconstraint::geq,
                                     dlpenalized ? -theta * deadlyloss : 0)); //c

        constraints->push_back(linearconstraint(
                                    dimupto(stage),
                                    linearconstraint::eq,
                                     stage ? fixedcost : 0)); //y

        constraints->push_back(linearconstraint(
                                    dimupto(stage),
                                    linearconstraint::eq,
                                     stage ? fixedcost : 0)); // z

        std::vector<double>& le = (constraints->rbegin()+3)->lhs;
        std::vector<double>& lc = (constraints->rbegin()+2)->lhs;
        std::vector<double>& ly = (constraints->rbegin()+1)->lhs;
        std::vector<double>& lz = constraints->rbegin()->lhs;

        lc[dimupto(stage)-cno] = 1;
        if(dlpenalized)
        {
            for(unsigned int s=0; s<=stage; s++)
                lc[dimupto(s)-zno] = -theta;
        }
        else
        {
            for(unsigned int s=0; s<=stage; s++)
                lc[dimupto(s)-zno] = -(laststage ? sigma : varrho);
        }
        le[dimupto(stage)-eno] = 1; // this stage e
        le[dimupto(stage)-sno] =-1; // this stage s

        lz[dimupto(stage)-zno] = 1.0; // z
        lz[dimupto(stage)-sno] = -xi[stage].p; // s
        if(stage)
            lz[dimupto(stage-1)-cno] = -1; // c

        if(stage)
        {
            assert(stageoffset(stage)-eno < dimupto(stage));
            le[stageoffset(stage)-eno]=-1; // former e

            unsigned int k = stageoffset(stage-1);

            for(unsigned int i=0; i<numsold;i++)
            {
                assert(k < dimupto(stage-1));
                le[k]=hR[i];
                lz[k]=xi[stage-1].m[i];
                k++;
            }
            for(int s = stage - 1, i=0; s>=0; s--, i++)
            {
                if(i<numfutures(s))
                {
                    assert(stageoffset(s)+numsold+i < dimupto(stage));

                    le[stageoffset(s)+numsold+i] = -1;
                    lz[stageoffset(s)+numsold+i] = -xi[s].q[i];
                }
            }
        }

        ly[dimupto(stage)-yno] = 1.0; // y
        ly[dimupto(stage)-sno] = -xi[stage].p; // s
        ly[dimupto(stage)-cno] = -1; // c

        if(!laststage)
        {
            unsigned int k = stageoffset(stage);

            for(unsigned int i=0; i<numsold;)
            {
                assert(k < dimupto(stage));
                ly[k++]=xi[stage].m[i++];
            }
            for(unsigned int i=0; i<numfutures(stage); i++)
                ly[k++] =  -xi[stage].q[i];
        }
   }
private:
    bool fz;
    bool flp;
    bool fzf;
    bool trivialproductionconstraint[numproducts];
};



class emitreesolutionstats: public treecallback
{
    enum foreachnodemode {eprob, estats};
public:

    emitreesolutionstats(const treesolution_ptr& sol) : fs(sol)
    {
    }


    double probdefault()
    {
        fimode = eprob;
        fi=0;
        fz.resize(fs->numstages());
        fp=0;
        fs->tp()->t()->foreachnode(this);
        return fp;
    }
private:

// state variables

    unsigned int fi;
    foreachnodemode fimode;
    std::vector<double> fz;
    prob fp;

    virtual void callback(const path& p)
    {
        if(fimode==eprob)
        {
            unsigned int s = p.size()-1;
            bool laststage = p.size()==fs->numstages();
            // tbd;
            for(unsigned int i=0; i<fs->sd(s); i++)
            {;
                //std::cout << "p.size()=" << p.size() << std::endl;
                //std::cout << "nx" << fs->numstages() << std::endl;

                //std::cout << "z=" << z << std::endl;
                int raoffset = s==0 || laststage ? 1 : 2;
                if(i==fs->sd(s)-emiproblem::zno-raoffset)
                {
                    fz[s]=fs->x(fi);
                    if(laststage)
                    {
                        double z = 0;
                        for(unsigned int j=0; j<fs->numstages();j++)
                            z += fz[j];

                        if(z > deadlyloss)
                            fp += fs->tp()->up(p);
                    }
               }
               fi++;
            }
        }
    }
private:
    treesolution_ptr fs;
};

class emitreemapping : public treemapping<emiatom>
{
public:
    emitreemapping(const treemapping_ptr<emixi>& xis, unsigned int pm) :
        treemapping<emiatom>(xis->t()), fxis(xis), fpm(pm)
      {}

    virtual scenario<emiatom> s(const path& p) const
    {
        double prevp=fpm*p0;
        std::vector<double> prevd(numsold);
        for(unsigned int i=0; i<numsold; i++)
            prevd[i]=Ed[i];

        scenario<emixi> exi = fxis->s(p);
        scenario<emiatom> r;

        for(unsigned int i=0; i<p.size(); i++)
        {
            r.push_back(emiatom(exi[i],prevp,prevd));
            prevp=r[r.size()-1].p;
            for(unsigned int j=0; j<numsold; j++)
            {
                prevd[j]=r[r.size()-1].d[j];
                assert(r.size()==5 || prevd[j]>0);
/*                if(minus)
                {
                    for(unsigned int i=0; i<numsold; i++)
                        std::cerr << d[i] << "(" << prevd[i] << "," << dr[i] << "),";
                    std::cerr << std::endl;
                    throw;
                }*/
            }
        }
        return r;
    }

    virtual emiatom operator()(const path& p) const
    {
        scenario<emiatom> s = this->s(p);
        return s[p.size()-1];
    }

private:
     treemapping_ptr<emixi> fxis;
     unsigned int fpm;
};

extern std::ostream& operator<<(std::ostream& os, const emiatom& e);


#endif // EMIMULTI_H
