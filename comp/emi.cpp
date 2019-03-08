#include <iostream>
#include <cmath>
#include "emimulti.h"
#include "tests.h"


std::ostream& operator<<(std::ostream& os, const emiatom& e)
{
    for(unsigned int i=0; ; i++)
    {
        os << e.d[i];
        os << ",";
        if(i==numsold-1)
            break;
    }
    for(unsigned int i=0; ; i++)
    {
        os << e.m[i];
        os << ",";
        if(i==numsold-1)
            break;
    }
    os << e.p << ",";
    for(unsigned int i=0; ; i++)
    {
        os << e.q[i];
        if(i==maxfutures-1)
            break;
        os << ",";
    }
    return os;
}






int main(int, char** /* int argc, char *argv[] */)
{
    const double c[] = {
        9.24,      // cuts
        92.7, // profiles
        20.6 // brams
    }; // contaminating parameter, the first actual ccordinate is computed elsewhere



    bool logsol = false; // if true then the solution is logged into a file
    bool pricedist2 = false; // two or three atom price distribution (three used in annor)


    const int T=4;
    bool z = true;
    bool dlpen = false; // deadly loss penalized

    double lora = whattodo < econtamination ? 0 : 0.2;
    double hira = whattodo < econtamination ? 1 : 0.2;
    double step = 0.1;

    double lomstdev = whattodo == emstac ? 0.05 : 0;
    double himstdev = whattodo == emstac ? 0.25 : 0;

    unsigned int lopm = whattodo < econtamination ? 0 : 1; // price multiplier for steps
    unsigned int hipm = whattodo < econtamination ? 8 : 1;

    double lonu = 0.08;
    double hinu = whattodo == econtamination ? 0.24 : 0.08;

    unsigned int lozf = 0;
    unsigned int hizf = whattodo < econtamination ? 1 : 0;

    std::vector<double> n3dist = {-n3atom, 0, n3atom};

    std::vector<double> n2dist = {-1,1};

    sys::set_log("mspp.log");

    std::ofstream os("results.csv");
    if(!os)
    {
        std::cerr << "Cannot create output file " << std::endl;
        throw;
    }

    if(whattodo==ednonstac)
    {
        alpha = alphanstac;
        beta = betanstac;
    }

    bool firsttime = true;
    for(double mstdev = lomstdev; mstdev < himstdev+0.001; mstdev+=0.1)
     for(double nu = lonu ; nu < hinu+0.001; nu+=0.08)
      for(unsigned int zerofutures = lozf; zerofutures <= hizf; zerofutures++)
        for(unsigned int pm=lopm; pm<hipm+1; pm= pm ? (pm*2) : 1)
        {
            std::cout << "Creating scenario tree for nu=" << nu <<
                         " mstdev=" << mstdev <<
                         " pm=" << pm << "..." << std::endl;

            sqV[4][4] = mstdev;

            scenariotree_ptr<emiatom> as;

            std::vector<std::vector<double>> first;
            first.push_back(xi0);

            std::vector<double>& udist = pricedist2 ? n2dist : n3dist;

            double mu = Exi3(udist,sqV[3][3]);

            std::vector<std::vector<double>> sd;
            if(whattodo==efixeddemand)
                sd = {{0}, {0}, {0}, udist, {0}};
            else if(whattodo==emstac )
                sd = {{-1,1}, {-1,1}, {-1,1}, udist,{-1,1}};
            else
                sd = {{-1,1}, {-1,1}, {-1,1}, udist, {0} };

            meanvardistribution middled({Exi1to3[0], Exi1to3[1],Exi1to3[2], mu,1} ,sqV,sd);


            std::vector<double> mlast = { 0,0,0,mu,1};

            meanvardistribution lastd(mlast,sqV, { {0}, {0},{0},udist,{0} } );
            std::vector<std::vector<double>> last = lastd.d();

            std::vector<std::vector<std::vector<double>>> dist;
            dist.push_back(first);

            if(whattodo == ebasic || whattodo == efixeddemand ||
                    whattodo == emstac || whattodo == ednonstac)
            {
                std::vector<std::vector<double>> middle = middled.d();

                for(unsigned int i=0; i<T-1; i++)
                    dist.push_back(middle);
                dist.push_back(last);
            }
            else if(whattodo == ewsscompute)
            {
                for(unsigned int i=0; i<T; i++)
                    dist.push_back(first);
            }
            else if(whattodo == econtamination)
            {
                std::vector<std::vector<double>> sd = {{-1,1}, {-1,1}, {-1,1}, udist,{0}};

                meanvardistribution middled({Exi1to3[0], Exi1to3[1],Exi1to3[2], mu,1} ,sqV,sd);
                std::vector<std::vector<double>> middle = middled.d();
                unsigned int ms = middle.size();
                unsigned int as = middle[0].size();

                middle.resize(2*ms);
                for(unsigned int i=0; i<ms; i++ )
                {
                    middle[i+ms] = middle[i];
                    unsigned int j=0;
                    for(;j<as-1;j++)
                        middle[i+ms][j]=middle[i][j]+c[j];
                    middle[i+ms][j]=middle[i][j];
                }

                for(unsigned int i=0; i<T-1; i++)
                    dist.push_back(middle);
                dist.push_back(last);

            }
            else assert(0);
            std::vector<std::vector<double>> ps;
            if(whattodo == econtamination)
            {
                ps.push_back({1});
                std::vector<double> prs(dist[1].size());
                unsigned int szhalf = prs.size() / 2;
                assert(szhalf == 24);

                for(unsigned int i=0; i<szhalf*2; i++)
                {
                    double coef = i < szhalf ? (1-nu) : nu;
                    prs[i]= 1.0 / ((double) szhalf) * coef;
                }


                for(unsigned int i=0; i<T-1; i++)
                    ps.push_back(prs);

                std::vector<double> lastprs(dist[T].size());
                for(unsigned int i=0; i< lastprs.size(); i++)
                    lastprs[i]= 1.0 / (double) lastprs.size();
                ps.push_back(lastprs);
            }

            idtreemapping_ptr<emixi> is(new idtreemapping<emixi>(dist));
            treemapping_ptr<emiatom> es(new emitreemapping(is,pm));


            treeprobability_ptr up(whattodo != econtamination ?
                                     (treeprobability*) new uniformtreeprobability(es->t())
                                    : (treeprobability*) new idtreeprobability(es->t(),ps));

            as.reset(new modularscenariotree<emiatom>(es, up));

//            scenariolister<emiatom> lister(as, os);
//            lister.list();
//            throw;


            linearproblem_ptr<emiatom> ep(new emiproblem(z,T, dlpen, zerofutures));

            for(double lambda=lora; lambda <= hira +0.000000001; lambda+=step)
            {
                std::cout << "p0=" << pm*p0 << ", lambda=" << lambda
                          << " futures=" << (zerofutures ? "no" : "yes") << "..." << std::endl;

                std::cout << "Creating meancvar variant of the problem ..." << std::endl;

                linearproblem_ptr<emiatom> cvp(new mmpcvarproblem<emiatom>(ep,0.05,lambda));

                if(firsttime)
                {
                    os << "p0,lambda,futures_excluded, optimal,";
                    for(unsigned int i=0;;)
                    {
                        os << cvp->varname(i);
                        if(++i==cvp->totaldim())
                            break;
                        os << ",";
                    }

                    os << ",pd" << std::endl;
                    firsttime = false;
                }
                os << pm*p0 << "," << lambda << ","
                   << (zerofutures ? 1 : 0) << ",";

                //csvlpsolver cps; // for debugging
                cplexlpsolver cps;

                std::cout << "Solving the problem..." << std::endl;

                biglpsolution<emiatom> b(cvp,as);

                treesolution_ptr s;
                double ov;

                b.solve(cps,s,ov);
                std::cout << "solved." << std::endl << std::endl;

                if(logsol)
                {
                    std::ostringstream solname;
                    solname << "sol_p0=" << pm*p0 << "_lambda=" << lambda << " .csv";

                    std::ofstream xls(solname.str());
                    xls << "path,";
                    for(unsigned int i=0;;)
                    {
                        xls << cvp->varname(i);
                        if(++i==cvp->totaldim())
                            break;
                        xls << ",";
                    }
                    xls << std::endl;
                    s->list(xls);
                }

                os << ov << ",";
                std::vector<double> E, var;
                s->stats(E,var);
                const double tol = 1e-7;
                for(unsigned int i=0; i<E.size(); i++)
                    E[i] = fabs(E[i]) < tol ? 0 : E[i];

                output(E,os);
                emitreesolutionstats ss(s);

                os << "," << ss.probdefault() << std::endl;
            }
        }
}


// tbd chytat chybu cplexu

