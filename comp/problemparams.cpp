#include "emimulti.h"

// paste from "problem params.2c++"

// (I-A)^{-1}
const double R[numproducts][numproducts]=
{		{	1	,	1.494	,	0	,	0	,	0	},
        {	0	,	1	,	0	,	0	,	0	},
        {	0	,	0	,	1	,	0	,	0	},
        {	0.645	,	0.96363	,	1.129	,	1	,	0	},
        {	0.48246	,	0.72079524	,	0.844492	,	0.748	,	1	}};

// margins EUR/1000t
const double	m	[numproducts]={	43350	,	230210	,	109430	,	19530	,	0	};

// production limits (in 1000t)
const double	w	[numproducts]={	750	,	30	,	120	,	900	,	680	};

// mean demand (1000t)
const double	Ed	[numproducts]={	510	,	28	,	90	,	20	,	0	};

// log conveniences
const double	Ey	[maxfutures]={	0.00974885	,	0.0194977	,	0.02924655	,	0.0389954			};

// emisions (t/1000t final)
const double	hR	[numproducts]={	258.941226	,	423.658191644	,	450.3870452	,	359.5988	,	113.1	};

// free alocated amounts (2017,18,19,20)
const double	r	[maxT]={	134847	,	132218	,	129555	,	126884	};

// spots
double p0 = 		6.54	;


// vector	(d_cuts,d_profiles,d_brams,r)

// mean
std::vector<double> Exi1to3=		{	28,	90,	20}	;
// variance
std::vector<std::vector<double>> sqV=
{		{	10.0465,	2.87088,	0,	0,	0	},
        {	2.87088,	14.1596,	0,	0,	0	},
        {	0,	0,	5.84715788292407,	0,	0	},
        {	0,	0,	0,	0.43854,	0	},
        {	0,	0,	0,	0,
0.15
            //0.55788
        }};

const double fixedcost =43282037;
const double deadlyloss = 10548816;


// atom of 3 point distribution matching N(0,1)
const double n3atom=		1.22474487139159	;

// end of paste

const std::vector<double> xi0 = { Ed[1],Ed[2],Ed[3], 0, 1 };

double alphanstac=	0.5	;
double betanstac=	1.8384	;


double alpha=	0;
double beta=	1;
