#define ltstep  	1d-10
#define smdist  	5d-12

#define dnc    		1.46d0
#define dcc    		1.51d0
#define dcn    		1.33d0

#define dcaca  		3.8d0
#define dtie   		2.41d0
#define dtie2  		2.45d0

#define sqz1   		1.1436d0
#define sqz2   		0.88d0
#define sqz3   		0.87829d0
#define sqz4   		0.8d0
#define sqz5   		0.7713d0

!#define sqz6                   0.7607d0
!sqz6_new based upon position
!#define sqz6            0.89552d0
!#define sqz7                   0.793d0
!sqz7_new based upon position
!#define sqz7            0.93333d0
!#define sqz8                   1.0956d0
!sqz8_new based upon position
!#define sqz8            1.375d0
!#define sqz9                   1.1244d0
!sqz9_new based upon position
!#define sqz9            1.3125d0
!#define sqz10                  0.9259d0
!sqz10_new based upon position
!#define sqz10           1.1619354d0

#define sqz11  		1.d0

#define n_b_hydro 	3
#define n_b_hbond 	3

#define del    		0.02375d0
#define rl_const 	1.2d0

#define numchains	(nop1/numbeads1)+(nop2/numbeads2)

#ifndef n_wrap
#define n_wrap 		1
#endif

#if n_wrap == 1
#define n_nab_cell      13
#else
#define n_nab_cell 	62
#endif

#define xrepuls1        40
#define xrepuls2        50

#define maxnbs 2000

#ifdef runs
#define	numsheets 3
#endif
