#include <iostream>
#include <stdio.h>
#include <string>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <chrono>
#include <string>
//#include <random>
#include "include/armadillo"
using namespace arma;
using namespace std;
//extern unsigned int seed;
//extern default_random_engine generator;
//uniform_int_distribution<unsigned int> unifage(5400,23400);
//uniform_real_distribution<double> unif01(0.0,1.0);
class Nodeman
{
public:
    unsigned int name;
    char gender[5] = "male";
    unsigned int age;
    vec status = zeros<vec>(3);       //double type
    umat periods = zeros<umat>(3, 4); //unsigned type uword
    vec resistance = zeros<vec>(3);   //double type
    unsigned int behaviour;
    unsigned int analrole;
    unsigned int anal_condom;
    unsigned int groupsex;
    double acqrate;
    unsigned int hiv = 0;
    unsigned int testday = 0;
    umat RhisPar = zeros<umat>(1, 4);
    umat ChisPar = zeros<umat>(1, 4);
    umat GhisPar = zeros<umat>(1, 4);
    Nodeman();
    Nodeman(unsigned int id, long long rseed);
    void Initializer(unsigned int id, long long rseed);
    static Nodeman *makeNodeman();
};
unsigned int getstopday(unsigned int t1, unsigned int t2, unsigned int t3);
void uvec_push(arma::uvec &v, unsigned int value);
void umat_push(arma::umat &v, unsigned int value1, unsigned int value2);
void rirj_partnership(Nodeman &mri, Nodeman &mrj, mat &eventsij, unsigned int &stopDayR, uword &ratri, uword &ratrj, unsigned int today, long long rseed);
void cicj_partnership(Nodeman &mci, Nodeman &mcj, mat &eventsij, unsigned int &stopDayC, uword &catci, uword &catcj, unsigned int today, long long rseed);
void gigj_partnership(Nodeman &mgi, Nodeman &mgj, mat &eventsij, rowvec &pgi_ac, rowvec &pgj_ac, unsigned int today, long long rseed);
mat sch_sex(mat events, long long rseed);
void transAB(Nodeman &A, Nodeman &B, mat events, unsigned typem, unsigned int today, vec p_act, urowvec &trans_sites, long long rseed);
//void outbreak_simulation(vec p_act, unsigned runid, unsigned my_years, unsigned site_1stAMR,
//                         vec &output, uvec &output_infected, umat &trans_tree,umat &outbreak_d, bool process, bool background,
//                         bool reconsc, bool regular, bool casual, double par_rcr, default_random_engine generator);
void outbreak_simulation(vec p_act, unsigned runid, unsigned my_years, unsigned site_1stAMR,
                         vec &output, uvec &output_infected, umat &trans_tree, umat &outbreak_d, bool process, bool background,
                         bool reconsc, bool regular, bool casual,
                         double par_rcr, double par_rr, double par_dignrate, double par_curerate,
                         default_random_engine generator);

//unsigned seed = chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator;

int main(int argc, char *argv[])
{

    clock_t startTime, endTime;
    startTime = clock();
    unsigned Runid = atoi(argv[1]);
    unsigned batchsize = 1; //500

    char temp1[100];
    //bool reconsc = false;
    //bool regular = false;
    //bool casual = false;
    bool reconsc = (*argv[2] == 'T');
    bool regular = (*argv[3] == 'T');
    bool casual = (*argv[4] == 'T');
    double par_rr = atof(argv[5]);
    double par_rcr = atof(argv[6]); //trace rate of casual partners

    sprintf(temp1, "./output_files_pes/OB_pes_recommended%d_regular%d_casual%d_%1.1f_%1.1f_%04d.txt", reconsc, regular, casual, par_rr, par_rcr, Runid); /////char temp[250]
    ofstream fp_output;
    fp_output.open(temp1, ios::trunc);
    sprintf(temp1, "./output_files_pes/OB_pes_spreadtree_recommended%d_regular%d_casual%d_%1.1f_%1.1f_%04d.txt", reconsc, regular, casual, par_rr, par_rcr, Runid); /////char temp[250]
    ofstream fp_output_spreadtree;
    fp_output_spreadtree.open(temp1, ios::trunc);
    sprintf(temp1, "./output_files_pes/OB_pes_outbreak_recommended%d_regular%d_casual%d_%1.1f_%1.1f_%04d.txt", reconsc, regular, casual, par_rr, par_rcr, Runid); /////char temp[250]
    ofstream fp_outbreak;
    fp_outbreak.open(temp1, ios::trunc);

    vec p_act = {0.0847642, 0.00940469, 0.469055, 0.0589935, 0.00491175, 0.694099, 0.0125556, 0.00422359};

    unsigned myyear = 5;      //20 year run after induction myyear+20*background+2
    unsigned site_1stAMR = 4; //atoi(argv[2]);

    double par_dignrate = 0.95;
    double par_curerate = 0.95;
    for (unsigned runid = batchsize * (Runid - 1); runid < batchsize * (Runid); runid++)
    {
        //runid+(Runid-1)*batchsize;
        vec output(2);
        uvec output_infected(1);
        umat trans_tree(1, 2);
        umat outbreak_d(2, 1);
        output.zeros();
        output_infected.zeros();
        trans_tree.zeros();

        default_random_engine generator(runid + 1);
        outbreak_d.resize(2, 1);
        outbreak_d(0, 0) = runid;
        outbreak_d(1, 0) = runid;
        outbreak_simulation(p_act, runid, myyear, site_1stAMR, output, output_infected, trans_tree, outbreak_d, false, true,
                            reconsc, regular, casual,
                            par_rcr, par_rr, par_dignrate, par_curerate,
                            generator);
        //bool reconsc, bool regular, bool casual,double par_rcr)
        //void outbreak_simulation(vec p_act, unsigned runid, unsigned my_years, unsigned site_1stAMR,
        //                 vec &output, uvec &output_infected, umat &trans_tree,umat &outbreak_d, bool process, bool background,
        //                 bool reconsc, bool regular, bool casual,
        //                 double par_rcr, double par_rr, double par_dignrate, double par_curerate,
        //                 default_random_engine generator)
        for (unsigned i = 0; i < output.n_elem; i++)
        {
            fp_output << (double)output(i) << " ";
        }
        for (unsigned i = 0; i < output_infected.n_elem; i++)
        {
            fp_output << (double)output_infected(i) << " ";
        }
        fp_output << std::endl;

        for (unsigned i = 0; i < trans_tree.n_cols; i++)
        {
            fp_output_spreadtree << (int)runid << " ";
            for (unsigned j = 0; j < trans_tree.n_rows; j++)
            {
                fp_output_spreadtree << (int)trans_tree(j, i) << " ";
            }
            fp_output_spreadtree << std::endl;
        }
        fp_output_spreadtree << std::endl;
        for (unsigned j = 0; j < 2; j++)
        {
            for (unsigned i = 0; i < outbreak_d.n_cols; i++)
            {
                fp_outbreak << (int)outbreak_d(j, i) << " ";
            }
            fp_outbreak << std::endl;
        }

        output.resize(2);
        output_infected.resize(1);
        trans_tree.resize(1, 2);
    }
    fp_output.close();
    fp_output_spreadtree.close();
    fp_outbreak.close();
    endTime = clock();
    std::cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
    return 0;
}

void outbreak_simulation(vec p_act, unsigned runid, unsigned my_years, unsigned site_1stAMR,
                         vec &output, uvec &output_infected, umat &trans_tree, umat &outbreak_d, bool process, bool background,
                         bool reconsc, bool regular, bool casual,
                         double par_rcr, double par_rr, double par_dignrate, double par_curerate,
                         default_random_engine generator) //recommended screening, casual par_rcr
{
    //set random seed with runid
    uniform_int_distribution<long long> rseed(1, 4294967296);
    uniform_real_distribution<double> unif01(0, 1.0);
    //uniform_real_distribution<double> gS(2.0, 4.0);
    gamma_distribution<double> gS(3.01, 0.8615);
    gamma_distribution<double> gAU(206.4, 0.4118);
    gamma_distribution<double> gO(201.6, 0.4167);
    gamma_distribution<double> gR(3702.857, 0.0972);
    normal_distribution<double> rnorm(0, 1);
    arma_rng::set_seed(rseed(generator));

    char temp1[100];
    char temp2[100];
    sprintf(temp1, "./output_files_pes/Prevalence%03d%02d.txt", runid, 1); /////char temp[250]
    ofstream fp_pre;
    sprintf(temp2, "./output_files_pes/Incidence%03d%02d.txt", runid, 1); /////char temp[250]
    ofstream fp_ind;
    //sprintf(temp3, "outputrandom%03d%02d.txt", runid, 1); /////char temp[250]
    //ofstream fp_rnum;
    //fp_rnum.open(temp3,ios::app);
    if (process)
    {
        fp_pre.open(temp1, ios::trunc);
        fp_ind.open(temp2, ios::trunc);
    }

    const int N = 10000;
    unsigned par_day = 5 * 360;
    unsigned back_years = 6 * (background == true);                 //
    unsigned max_day = 360 * my_years + par_day + back_years * 360; //2 years run for network, 10 years run for status
    bool gs_trace = false;
    vec culture_rate(3);
    culture_rate(0) = 1; //use rate data
    culture_rate(1) = 1;
    culture_rate(2) = 1;

    bool Debuge = false;
    Nodeman *men = new Nodeman[N];

    mat prob_ac = zeros<mat>(N, 2);
    uvec scrday(N);
    scrday.zeros();
    uvec Gat(N);
    Gat.zeros();
    double U;

    //Initialization of agents
    uvec avR(N);
    avR.ones();
    uvec avC(N);
    avC.ones();
    uvec avG(N);
    avG.fill(100); //people who have group sex, the value is 0/1, and 100 means no group sex

    vec acqrate(N);
    acqrate.zeros();

    for (unsigned int i = 0; i < N; i++)
    {
        men[i].Initializer(i, rseed(generator));
        prob_ac(i, 0) = (double)men[i].acqrate / 180;
        prob_ac(i, 1) = 0;
        if (men[i].acqrate > 5)
        {
            men[i].groupsex = (double)23 / 55 * (men[i].acqrate - 5) + 1; // most 24 times
            U = unif01(generator);
            Gat(i) = ceil(-log(U) * (double)360 / men[i].groupsex) + 1; //get available for group sex
            avG(i) = 0;
            //fp_rnum<<"Group sex"<<U<<" "<<Gat(i)<< endl;
        }
        acqrate(i) = men[i].acqrate;
        U = unif01(generator);
        scrday(i) = ceil(U * 360); //ceil((randn() / 10 + 1) * men[i].testday);//need to be normal
        if (men[i].behaviour == 1)
        {
            avC(i) = 100;
        }
        else if (men[i].behaviour == 2)
        {
            avR(i) = 100;
        }
    } //End initialization of agents

    //fp_rnum<<men[10].acqrate<<' '<<men[100].acqrate<<' '<<men[100].acqrate<<endl;
    //fp_rnum<< "s0"<<unif01(generator)<<" "<<randu<double>()<<endl;

    uvec hrskid;
    uvec lrskid;
    hrskid = find(acqrate > 9);
    lrskid = find(acqrate <= 9);
    std::cout << "Begin program" << std::endl;
    uvec Rat(N);
    uvec Cat(N);
    uvec Ra = find(avR == 1);
    uvec Ca = find(avC == 1);

    uvec Patstypes;
    umat Patsid(0, 2);
    mat EvM0(0, 5);
    mat EvM1(0, 5);
    mat events(5, 2);
    uvec StopDay;
    StopDay.zeros();

    unsigned int ri, rj;
    unsigned int ci, cj;
    unsigned int m1, m2;
    unsigned int stopday;
    while (Ra.n_elem > 1 && Ca.n_elem > 1)
    {
        U = unif01(generator);
        if (U < 0.5) //Ra
        {
            //cout << "RA" << endl;
            U = unif01(generator);
            m1 = floor(Ra.n_elem * U);
            ri = Ra(m1);
            Ra.shed_row(m1);
            U = unif01(generator);
            m2 = floor(Ra.n_elem * U);
            rj = Ra(m2);
            Ra.shed_row(m2);
            if (unif01(generator) < 0.5)
                continue;
            //toatal partnership is zero, then insert at zeroth row.
            events.zeros();
            rirj_partnership(men[ri], men[rj], events, stopday, Rat(ri), Rat(rj), 0, rseed(generator));
            uvec_push(Patstypes, 1);
            umat_push(Patsid, ri, rj);
            EvM0.insert_rows(EvM0.n_rows, events.col(0).t());
            EvM1.insert_rows(EvM1.n_rows, events.col(1).t());
            uvec_push(StopDay, stopday);
            avR(ri) = 0;
            avR(rj) = 0;
        }
        else // Ca
        {
            //cout << "CA" << endl;
            U = unif01(generator);
            m1 = floor(Ca.size() * U);
            ci = Ca(m1);
            Ca.shed_row(m1);
            U = unif01(generator);
            m2 = floor(Ca.size() * U);
            cj = Ca(m2);
            Ca.shed_row(m2);
            if (unif01(generator) < 0.5)
                continue;
            events.zeros();
            cicj_partnership(men[ci], men[cj], events, stopday, Cat(ci), Cat(cj), 0, rseed(generator));

            uvec_push(Patstypes, 2);
            umat_push(Patsid, ci, cj);
            EvM0.insert_rows(EvM0.n_rows, events.col(0).t());
            EvM1.insert_rows(EvM1.n_rows, events.col(1).t());
            uvec_push(StopDay, stopday);

            avC(ci) = 0;
            avC(cj) = 0;
        };
    }
    //fp_rnum<< "s1"<<unif01(generator)<<" "<<randu<double>()<<endl;

    umat diagnosis(N, 4);
    diagnosis.zeros();
    umat laststatus(N, 4);
    laststatus.zeros();
    umat lastAMR(N, 4);
    lastAMR.zeros();
    umat diagnosis_AMR(N, 4);
    diagnosis_AMR.zeros();

    //To calculate Denton's incidence
    mat incd_num_test(N, 4);
    mat incd_lasttest(N, 4);
    incd_num_test.zeros();
    incd_lasttest.zeros();
    uvec id_test(N);
    id_test.zeros();
    rowvec ind_high(4);
    ind_high.zeros();
    uvec last_testday(N);
    last_testday.zeros(N);

    //To calculate prevalence
    urowvec all_pre(4);
    urowvec h_pre(4);
    urowvec l_pre(4);
    //To calculate prevalence of AMR
    urowvec all_AMR(4);
    urowvec h_AMR(4);
    urowvec l_AMR(4);
    //To calculate real incidence
    urowvec all_ind(4);
    urowvec h_ind(4);
    urowvec l_ind(4);
    bool site_treated;
    urowvec site_treated_all(3);
    double detectAMR_prob = 0;
    bool detectAMR = false;
    bool detected1st_test = false;
    bool detected1st_treat = false;

    urowvec cur_stats(4);
    urowvec cur_AMR(4);
    urowvec ind_temp(4);
    urowvec AMR_temp(4);
    urowvec total_AMR(4);
    total_AMR.zeros(); // acumulated AMR cases
    uvec outbreak_d_temp(2);
    vec info_1stinfection(7);
    info_1stinfection.zeros();
    vec info_1stdetected(19);
    info_1stdetected.zeros();

    unsigned int rp; //new enter person
    unsigned int NumPar6m = 0;
    uvec trace_partners(1);
    vec Num_tt(2);
    Num_tt.zeros();

    urowvec total_trans_sites(8);
    total_trans_sites.zeros();
    urowvec trans_sites(8);
    std::cout << "Burn-in period";
    for (unsigned int d = 1; d < max_day + 1; d++)
    {
        if (d >= par_day + back_years * 360 && (d - par_day - back_years * 360 - 1) % 360 == 0)
        {
            std::cout << " " << std::endl;
            std::cout << "Year=" << (d - par_day - back_years * 360 - 1) / 360 + 1 << endl;
        }

        //execute sex acts and perform transimission
        //uvec IJ = find(TypeM > 0);
        uvec ij;
        vec ij0_status(3);
        vec ij1_status(3);
        for (unsigned int n = 0; n < Patsid.n_rows; n++)
        {
            ij = Patsid.row(n).t(); //ind2sub(arma::size(TypeM), IJ(n));
            events.zeros();
            EvM0.row(n).t();
            events.col(0) = EvM0.row(n).t();
            events.col(1) = EvM1.row(n).t();
            //cout << "Sex acts: Begin" << endl;
            //individual has symptomatic NG, which is uncomfortable for sexual acts
            if ((men[ij(0)].status(1) == 1 && d > men[ij(0)].periods(1, 1)) || (men[ij(0)].status(2) == 1 && d > men[ij(0)].periods(2, 1)) ||
                (men[ij(1)].status(1) == 1 && d > men[ij(1)].periods(1, 1)) || (men[ij(1)].status(2) == 1 && d > men[ij(1)].periods(2, 1)))
            {
                //cout << "Continue other couples" << endl;
                continue;
            }
            // both individuals are good so schedule next day events.
            if (all(men[ij(0)].status <= 0) && all(men[ij(1)].status <= 0))
            {
                events = sch_sex(events, rseed(generator)); //Schedule sex for next day
                EvM0.row(n) = events.col(0).t();
                EvM1.row(n) = events.col(1).t();
            }
            else
            { //TranAB
                //cout << "TransAB: Begin" << endl;
                ij0_status = men[ij(0)].status;
                ij1_status = men[ij(1)].status;
                trans_sites.zeros();
                transAB(men[ij(0)], men[ij(1)], events, Patstypes(n), d, p_act, trans_sites, rseed(generator));
                total_trans_sites = total_trans_sites + trans_sites;
                //if (!detected1st_treat && !detected1st_test && d > par_day + back_years * 360)
                if (d < par_day + back_years * 360 + 360 && d > par_day + back_years * 360)
                {
                    if (all(ij0_status == 0) && any(men[ij(0)].status > 0))
                    {
                        umat_push(trans_tree, ij(1), ij(0));
                    }
                    else if (all(ij1_status == 0) && any(men[ij(1)].status > 0))
                    {
                        umat_push(trans_tree, ij(0), ij(1));
                    }
                }

                events = sch_sex(events, rseed(generator)); //Schedule sex for next day
                EvM0.row(n) = events.col(0).t();
                EvM1.row(n) = events.col(1).t();
            }
        }

        //if (true)
        //{
        //    std::cout << "Recovery and treatment: Begin" << std::endl;
        //}
        //recovery and treatment d
        //assume that three anatomic sites are independatn.
        for (unsigned int i = 0; i < N; i++)
        {
            //For sympotmatical urethral and rectual, treatment is a detection for resistance as it will for sure have treatment failure.

            //---------------
            detectAMR = false; //indicator of whether response
            site_treated_all.zeros();
            if ((men[i].status(1) == 1 && men[i].periods(1, 2) == d) || (men[i].status(2) == 1 && men[i].periods(2, 2) == d))
            {
                last_testday(i) = d;
                site_treated_all = (men[i].resistance.t() == 0) * 1;                        //(O,U,R) resistant or not, not resistant, treat (1,1,1)
                scrday(i) = scrday(i) + ceil((rnorm(generator) / 10 + 1) * men[i].testday); //put backward next test

                if (all(men[i].resistance == 0))
                {
                    site_treated_all.ones();
                }
                else if ((men[i].status(1) == 1 && men[i].resistance(1) == 1) || (men[i].status(2) == 1 && men[i].resistance(2) == 1)) // symptomatic site
                {
                    site_treated_all.ones();
                    detectAMR = true;
                }
                else // other sites
                {
                    detectAMR_prob = 1 - pow(1 - culture_rate(0), men[i].resistance(0)) * pow(1 - culture_rate(1), men[i].resistance(1)) * pow(1 - culture_rate(2), men[i].resistance(2));
                    detectAMR = unif01(generator) < detectAMR_prob;
                    if (detectAMR)
                    {
                        site_treated_all.ones();
                    }
                }

            } //change
            if (!detected1st_treat && !detected1st_test && detectAMR)
            {
                detected1st_treat = true;
                info_1stdetected(0) = men[i].behaviour;
                info_1stdetected(1) = men[i].acqrate;
                info_1stdetected(2) = men[i].testday;
                info_1stdetected(3) = men[i].status(0);
                info_1stdetected(4) = men[i].status(1);
                info_1stdetected(5) = men[i].status(2);
                info_1stdetected(6) = men[i].resistance(0);
                info_1stdetected(7) = men[i].resistance(1);
                info_1stdetected(8) = men[i].resistance(2);
                info_1stdetected(9) = d - par_day - back_years * 360;
                all_AMR.zeros();
                //total_AMR.zeros();
                for (unsigned n = 0; n < N; n++)
                {
                    AMR_temp.zeros();
                    cur_AMR(0) = (men[n].resistance(0) > 0);
                    cur_AMR(1) = (men[n].resistance(1) > 0);
                    cur_AMR(2) = (men[n].resistance(2) > 0);
                    cur_AMR(3) = any(men[n].resistance > 0);
                    all_AMR = all_AMR + cur_AMR;
                    AMR_temp = cur_AMR > lastAMR.row(n);
                    diagnosis_AMR.row(n) = diagnosis_AMR.row(n) + AMR_temp;
                    //total_AMR = total_AMR + AMR_temp;
                }
                for (unsigned s = 0; s < 4; s++)
                {
                    info_1stdetected(10 + s) = all_AMR(s);
                    info_1stdetected(14 + s) = total_AMR(s);
                }
                if (men[i].RhisPar.n_rows > 1)
                {
                    trace_partners.resize(1);
                    trace_partners(0) = men[i].RhisPar(men[i].RhisPar.n_rows - 1, 1); //current regular parnter
                    info_1stdetected(18) = any(men[trace_partners(0)].resistance > 0);
                }
                else
                {
                    info_1stdetected(18) = 0;
                }

                //Implementation of recommened screening
                if (reconsc)
                {
                    for (unsigned n = 0; n < N; n++)
                    {
                        if (men[n].acqrate > 9)
                        {
                            men[n].testday = 90;
                        }
                        else
                        {
                            if (men[n].testday > 360)
                            {
                                men[n].testday = 360;
                            }
                        }
                    }
                }

            } // the first cases,the responese starts

            if (any(site_treated_all > 0))
            {
                //id_test = 0 not be tested; =1, had a negative test; =2 had a postive test;
                incd_lasttest(i, 0) = men[i].status(0) > 0;
                incd_lasttest(i, 1) = men[i].status(1) > 0;
                incd_lasttest(i, 2) = men[i].status(2) > 0;
                incd_lasttest(i, 3) = any(men[i].status > 0);
                if (id_test(i) == 1)
                {
                    incd_num_test.row(i) = incd_num_test.row(i) + incd_lasttest.row(i);
                }
                id_test(i) = 2;
                //}

                //if any sites need to be treated
                //treatment fails and success
                if (unif01(generator) < par_curerate) //cured
                {
                    for (unsigned int s = 0; s < 3; s++)
                    {
                        if (site_treated_all(s) == 1 && men[i].status(s) > 0) //&& men[i].resistance(s) == 0) //whether the screending day should change
                        {
                            men[i].periods(s, 0) = 0;
                            men[i].periods(s, 1) = 0;
                            men[i].periods(s, 2) = d;
                            men[i].periods(s, 3) = d + 7;
                        }
                    }
                }
                else //fails symptoms will keep on for one weeks, but no sexual contact with others.
                {
                    if (men[i].status(1) == 1 && site_treated_all(1) == 1)
                    {
                        men[i].periods(1, 0) = 0;
                        men[i].periods(1, 1) = d + 7;
                        men[i].periods(1, 2) = d + 14;
                    }

                    if (men[i].status(2) == 1 && site_treated_all(2) == 1)
                    {
                        men[i].periods(2, 0) = 0;
                        men[i].periods(2, 1) = d + 7;
                        men[i].periods(2, 2) = d + 14;
                    }
                }
                Num_tt(0) = Num_tt(0) + 1;
            }

            ////}////////confirm
            // Implemention of contact tracing of AMR cases
            if (regular)
            {
                trace_partners.resize(1);
                if (detectAMR && men[i].behaviour != 2 && men[i].RhisPar.n_rows > 1)
                {
                    trace_partners(0) = men[i].RhisPar(men[i].RhisPar.n_rows - 1, 1); //current regular parnter
                    if (unif01(generator) < par_rr && last_testday(trace_partners(0)) + 14 < d)
                    {
                        scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //come to test
                        //cout<<scrday(trace_partners(0))<<endl;
                    }
                }
            }
            if (casual)
            {
                if (detectAMR && men[i].behaviour != 1 && men[i].ChisPar.n_rows > 1)
                {
                    for (unsigned k = 0; k < men[i].ChisPar.n_rows - 1; k++)
                    {
                        if (men[i].ChisPar(men[i].ChisPar.n_rows - 1 - k, 3) > d - 60)
                        {
                            trace_partners(0) = men[i].ChisPar(men[i].ChisPar.n_rows - 1 - k, 1);
                            if (unif01(generator) < par_rcr && last_testday(trace_partners(0)) + 14 < d)
                            {
                                scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //come to test
                            }
                        }
                    }
                }
                if (detectAMR && men[i].behaviour != 1 && men[i].GhisPar.n_rows > 1 && gs_trace)
                {
                    for (unsigned k = 0; k < men[i].GhisPar.n_rows - 1; k++)
                    {
                        if (men[i].GhisPar(men[i].GhisPar.n_rows - 1 - k, 3) > d - 60)
                        {
                            trace_partners(0) = men[i].GhisPar(men[i].GhisPar.n_rows - 1 - k, 1);
                            if (unif01(generator) < par_rcr && last_testday(trace_partners(0)) + 14 < d)
                            {
                                scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //come to test
                            }
                        }
                    }
                }
            }
            ////}/////
            //recovery
            for (unsigned int s = 0; s < 3; s++)
            {

                if (men[i].periods(s, 2) == d)
                {
                    men[i].status(s) = -1;
                    men[i].resistance(s) = 0;
                    men[i].periods(s, 1) = 0;
                    men[i].periods(s, 2) = 0;
                }
                else if (men[i].periods(s, 3) == d)
                {
                    men[i].status(s) = 0;
                    men[i].periods(s, 3) = 0;
                }
            }

            //-------------
        } //end for i=1...N
        //if (true)
        //{
        //    std::cout << "Recovery and treatment: End" << std::endl;
        //    std::cout << "Separation on d and formation: Begin" << std::endl;
        //}
        //Separation on d and formation for d+1
        //Regular partnership separantion

        uvec RIJ = find(StopDay == d && Patstypes == 1);
        uvec rij;
        if (!RIJ.is_empty())
        {
            for (unsigned i = 0; i < RIJ.n_elem; i++)
            {
                rij = Patsid.row(RIJ(i) - i).t();
                Patsid.shed_row(RIJ(i) - i);
                StopDay.shed_row(RIJ(i) - i);
                Patstypes.shed_row(RIJ(i) - i);
                EvM0.shed_row(RIJ(i) - i);
                EvM1.shed_row(RIJ(i) - i);

                Rat(rij(0)) = 0;
                Rat(rij(1)) = 0;
                avR(rij(0)) = 1;
                avR(rij(1)) = 1;
            }
        }
        //Casual partnership separantion
        uvec CIJ = find(StopDay == d && Patstypes == 2);
        uvec cij;
        if (!CIJ.is_empty())
        {
            for (unsigned i = 0; i < CIJ.n_elem; i++)
            {

                cij = Patsid.row(CIJ(i) - i).t();
                Patsid.shed_row(CIJ(i) - i);
                StopDay.shed_row(CIJ(i) - i);
                Patstypes.shed_row(CIJ(i) - i);
                EvM0.shed_row(CIJ(i) - i);
                EvM1.shed_row(CIJ(i) - i);
                avC(cij(0)) = 2;
                avC(cij(1)) = 2;
                prob_ac(cij(0), 1) = d;
                prob_ac(cij(1), 1) = d;
            }
        }
        uvec CI = find(Cat == d);
        if (!CI.is_empty())
        {
            for (unsigned i = 0; i < CI.n_elem; i++)
            {
                avC(CI(i)) = 1;
            }
        }

        uvec GIJ = find(Patstypes == 3);
        uvec gij;
        if (!GIJ.is_empty())
        {
            for (unsigned i = 0; i < GIJ.n_elem; i++)
            {
                gij = Patsid.row(GIJ(i) - i).t();
                Patsid.shed_row(GIJ(i) - i);
                StopDay.shed_row(GIJ(i) - i);
                Patstypes.shed_row(GIJ(i) - i);
                EvM0.shed_row(GIJ(i) - i);
                EvM1.shed_row(GIJ(i) - i);
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "Separation on d and formation: End" << std::endl;
        //    std::cout << "Replacement and aging: Begin" << std::endl;
        //}
        //replacement and aging d+1, and probability update.
        for (unsigned n = 0; n < N; n++)
        {

            if (men[n].age + 1 > 23400)
            {                                                          //65 years od 23400
                uvec pri = find(Patsid.col(0) == n && Patstypes == 1); //find(Patsid.col(1)==n && Patstypes==1));
                uvec ri;

                if (!pri.is_empty())
                {
                    ri = Patsid(pri(0), 1);
                }
                else
                {
                    pri = find(Patsid.col(1) == n && Patstypes == 1);

                    if (!pri.is_empty())
                        ri = Patsid(pri(0), 0);
                }
                //ri.print();
                if (!ri.is_empty())
                {
                    Patsid.shed_row(pri(0));
                    StopDay.shed_row(pri(0));
                    Patstypes.shed_row(pri(0));
                    EvM0.shed_row(pri(0));
                    EvM1.shed_row(pri(0));
                    Rat(ri(0)) = 0;
                    avR(ri(0)) = 1;
                    men[ri(0)].RhisPar(men[ri(0)].RhisPar.n_rows - 1, 3) = d;
                }
                uvec pci = find(Patsid.col(0) == n && Patstypes == 2);
                uvec ci;
                if (!pci.is_empty())
                {
                    ci = Patsid(pci(0), 1);
                }
                else
                {
                    pci = find(Patsid.col(1) == n && Patstypes == 2);
                    if (!pci.is_empty())
                        ci = Patsid(pci(0), 0);
                }
                if (!ci.is_empty())
                {
                    Patsid.shed_row(pci(0));
                    StopDay.shed_row(pci(0));
                    Patstypes.shed_row(pci(0));
                    EvM0.shed_row(pci(0));
                    EvM1.shed_row(pci(0));

                    men[ci(0)].ChisPar(men[ci(0)].ChisPar.n_rows - 1, 3) = d;
                    avC(ci(0)) = 2;
                }
                men[n].age = 5760;
                //men[n].order_id = N + rp;
                rp++;
                men[n].status.zeros();
                men[n].periods.zeros();
                men[n].resistance.zeros();
                men[n].RhisPar = zeros<umat>(1, 4);
                men[n].ChisPar = zeros<umat>(1, 4);
                men[n].GhisPar = zeros<umat>(1, 4);
                avR(n) = (men[n].behaviour == 1) * 1 + (men[n].behaviour == 3) * 1 + (men[n].behaviour == 2) * 2;
                avC(n) = (men[n].behaviour == 1) * 2 + (men[n].behaviour == 3) * 1 + (men[n].behaviour == 2) * 1;
                avG(n) = (avG(n) < 10) * 1;
                Rat(n) = 0;
                Cat(n) = 0;
                Gat(n) = 0;
                prob_ac(n, 0) = (double)men[n].acqrate / 180;
                prob_ac(n, 1) = 0;
            }
            else
            {
                men[n].age = men[n].age + 1;
            }
            uvec past_parid;
            past_parid = find(men[n].RhisPar.col(3) <= d - 180);
            if (past_parid.n_elem > 1 && d > 180)
            {
                men[n].RhisPar.shed_rows(1, past_parid(past_parid.n_elem - 1));
            }
            past_parid = find(men[n].ChisPar.col(3) <= d - 180);

            if (past_parid.n_elem > 1 && d > 180)
            {
                men[n].ChisPar.shed_rows(1, past_parid(past_parid.n_elem - 1));
            }

            past_parid = find(men[n].GhisPar.col(3) <= d - 180);

            if (past_parid.n_elem > 1 && d > 180)
            {
                men[n].GhisPar.shed_rows(1, past_parid(past_parid.n_elem - 1));
            }
            NumPar6m = men[n].ChisPar.n_rows + men[n].GhisPar.n_rows - 2;
            if (NumPar6m > men[n].acqrate)
            {
                prob_ac(n, 0) = 0;
            }
            else
            {
                unsigned dayC = 0, dayG = 0;
                if (men[n].ChisPar.n_rows > 1)
                {
                    dayC = men[n].ChisPar(1, 3);
                }
                if (men[n].GhisPar.n_rows > 1)
                {
                    dayG = men[n].GhisPar(1, 3);
                }
                if (dayC > 0 && dayG == 0)
                {
                    prob_ac(n, 0) = (double)(men[n].acqrate - NumPar6m) / (180 - d + dayC);
                }
                else if (dayC == 0 && dayG > 0)
                {
                    prob_ac(n, 0) = (double)(men[n].acqrate - NumPar6m) / (180 - d + dayG);
                }
                else if (dayC > 0 && dayG > 0)
                {
                    prob_ac(n, 0) = (double)(men[n].acqrate - NumPar6m) / (180 - d + min({dayC, dayG}));
                }
                else
                {
                    prob_ac(n, 0) = (double)1 / (180 / men[n].acqrate - d + prob_ac(n, 1));
                    prob_ac(n, 0) = (double)1 * (prob_ac(n, 0) < 0) + prob_ac(n, 0) * (prob_ac(n, 0) > 0);
                    //why here is 1?
                }
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "Replacement and aging: End" << std::endl;
        //    std::cout << "Formation  of partnership: Begin" << std::endl;
        //}
        //Formation of group sex
        uvec GI = find(Gat == d);
        if (!GI.is_empty())
        {
            for (unsigned i = 0; i < GI.n_elem; i++)
            {
                avG(GI(i)) = 1;
            }
            //avG.rows(GI)=1;
        }
        //cout <<"Begin group size 1" << endl;
        uvec Ga = find(avG == 1 && prob_ac.col(0) > 0);
        mat ProbGa = prob_ac.rows(Ga);

        uword gi, gj;
        unsigned lenGa01 = ceil(Ga.n_elem * 0.10);

        while (Ga.n_elem > lenGa01)
        {
            //cout <<"Begin group size 2" << endl;
            unsigned GroupSize = floor(3 - 1.9 * log(unif01(generator)));
            vector<uword> Group;
            if (Ga.n_elem < GroupSize)
            {
                break;
            }
            else
            {
                for (unsigned i = 0; i < GroupSize; i++)
                {
                    //cout <<"Begin group size" << endl;
                    m1 = min(find(unif01(generator) < cumsum(ProbGa.col(0) / sum(ProbGa.col(0)))));
                    Group.push_back(Ga(m1));
                    Gat(Ga(m1)) = d + ceil(-log(unif01(generator)) * (double)360 / men[Ga(m1)].groupsex) + 1;
                    avG(Ga(m1)) = 0;
                    Ga.shed_row(m1);
                    ProbGa.shed_row(m1);
                }
                for (unsigned i = 0; i < GroupSize - 1; i++)
                {
                    for (unsigned j = i + 1; j < GroupSize; j++)
                    {
                        gi = Group[i];
                        gj = Group[j];
                        if (unif01(generator) < 0.8)
                        {
                            rowvec prob_ac_gi = prob_ac.row(gi);
                            rowvec prob_ac_gj = prob_ac.row(gj);
                            events.zeros();
                            gigj_partnership(men[gi], men[gj], events, prob_ac_gi, prob_ac_gj, d, rseed(generator));
                            uvec_push(Patstypes, 3);
                            umat_push(Patsid, gi, gj);
                            EvM0.insert_rows(EvM0.n_rows, events.col(0).t());
                            EvM1.insert_rows(EvM1.n_rows, events.col(1).t());
                            uvec_push(StopDay, d + 1);

                            prob_ac.row(gi) = prob_ac_gi;
                            prob_ac.row(gj) = prob_ac_gj;
                        }
                    }
                }
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "End formation  of group sex partnership" << std::endl;
        //}
        Ra = find(avR == 1);
        Ca = find(avC == 1 && prob_ac.col(0) > 0);

        while (Ra.n_elem > 1 || Ca.n_elem > 1)
        {

            double prob_re = (double)Ra.n_elem / (Ra.n_elem + Ca.n_elem);
            U = unif01(generator);
            //cout<<U<<endl;
            //cout<<prob_re<<endl;
            if (Ra.n_elem > 1 && U < prob_re)
            {
                U = unif01(generator);
                m1 = floor(Ra.n_elem * U);
                ri = Ra(m1);
                Ra.shed_row(m1);
                U = unif01(generator);
                m2 = floor(Ra.n_elem * U);
                rj = Ra(m2);
                Ra.shed_row(m2);

                if (unif01(generator) < 0.5)
                    continue;
                // write below into function rirj_partnership
                //rirj_partnership(men[ri], men[rj], EvM(ri, rj), stopday, Rat(ri), Rat(rj), d);

                events.zeros();
                rirj_partnership(men[ri], men[rj], events, stopday, Rat(ri), Rat(rj), d, rseed(generator));
                uvec_push(Patstypes, 1);
                umat_push(Patsid, ri, rj);
                EvM0.insert_rows(EvM0.n_rows, events.col(0).t());
                EvM1.insert_rows(EvM1.n_rows, events.col(1).t());
                uvec_push(StopDay, stopday);

                avR(ri) = 0;
                avR(rj) = 0;
            }
            else if (Ca.n_elem > 1 && U > prob_re)
            {
                U = unif01(generator);
                m1 = floor(Ca.n_elem * U);
                ci = Ca(m1);
                Ca.shed_row(m1);
                U = unif01(generator);
                m2 = floor(Ca.n_elem * U);
                cj = Ca(m2);
                Ca.shed_row(m2);
                if (unif01(generator) < 0.5)
                    continue;
                // write below into function rirj_partnership
                events.zeros();
                cicj_partnership(men[ci], men[cj], events, stopday, Cat(ci), Cat(cj), d, rseed(generator));

                uvec_push(Patstypes, 2);
                umat_push(Patsid, ci, cj);
                EvM0.insert_rows(EvM0.n_rows, events.col(0).t());
                EvM1.insert_rows(EvM1.n_rows, events.col(1).t());
                uvec_push(StopDay, stopday);

                avC(ci) = 0;
                avC(cj) = 0;
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "Formation of partnership: End" << std::endl;
        //number of cases
        //calculation of indicence  and prevalence
        //    std::cout << "Calculation of indicence  and prevalence: Begin" << std::endl;
        //}
        h_pre.zeros();
        l_pre.zeros();
        all_pre.zeros();
        h_AMR.zeros();
        l_AMR.zeros();
        all_AMR.zeros();
        if (d > par_day)
        {
            for (unsigned i = 0; i < hrskid.n_elem; i++)
            {
                ind_temp.zeros();
                cur_stats(0) = (men[hrskid(i)].status(0) > 0);
                cur_stats(1) = (men[hrskid(i)].status(1) > 0);
                cur_stats(2) = (men[hrskid(i)].status(2) > 0);
                cur_stats(3) = any(men[hrskid(i)].status > 0);
                h_pre = h_pre + cur_stats;
                all_pre = all_pre + cur_stats;
                ind_temp = (cur_stats > laststatus.row(hrskid(i)));
                diagnosis.row(hrskid(i)) = diagnosis.row(hrskid(i)) + ind_temp;

                AMR_temp.zeros();
                cur_AMR(0) = (men[hrskid(i)].resistance(0) > 0);
                cur_AMR(1) = (men[hrskid(i)].resistance(1) > 0);
                cur_AMR(2) = (men[hrskid(i)].resistance(2) > 0);
                cur_AMR(3) = any(men[hrskid(i)].resistance > 0);
                h_AMR = h_AMR + cur_AMR;
                all_AMR = all_AMR + cur_AMR;
                AMR_temp = cur_AMR > lastAMR.row(hrskid(i));
                diagnosis_AMR.row(hrskid(i)) = diagnosis_AMR.row(hrskid(i)) + AMR_temp;

                total_AMR = total_AMR + AMR_temp;
            }
            for (unsigned i = 0; i < lrskid.n_elem; i++)
            {
                ind_temp.zeros();
                cur_stats(0) = (men[lrskid(i)].status(0) > 0);
                cur_stats(1) = (men[lrskid(i)].status(1) > 0);
                cur_stats(2) = (men[lrskid(i)].status(2) > 0);
                cur_stats(3) = any(men[lrskid(i)].status > 0);
                l_pre = l_pre + cur_stats;
                all_pre = all_pre + cur_stats;
                ind_temp = (cur_stats > laststatus.row(lrskid(i)));
                diagnosis.row(lrskid(i)) = diagnosis.row(lrskid(i)) + ind_temp;

                AMR_temp.zeros();
                cur_AMR(0) = (men[lrskid(i)].resistance(0) > 0);
                cur_AMR(1) = (men[lrskid(i)].resistance(1) > 0);
                cur_AMR(2) = (men[lrskid(i)].resistance(2) > 0);
                cur_AMR(3) = any(men[lrskid(i)].resistance > 0);
                l_AMR = l_AMR + cur_AMR;
                all_AMR = all_AMR + cur_AMR;
                AMR_temp = cur_AMR > lastAMR.row(lrskid(i));
                diagnosis_AMR.row(lrskid(i)) = diagnosis_AMR.row(lrskid(i)) + AMR_temp;
                total_AMR = total_AMR + AMR_temp;
            }

            if (fp_pre.is_open())
            {
                fp_pre << (int)d << " " << (double)h_pre(0) / hrskid.n_elem * 100 << " " << (double)h_pre(1) / hrskid.n_elem * 100 << " "
                       << (double)h_pre(2) / hrskid.n_elem * 100 << " " << (double)h_pre(3) / hrskid.n_elem * 100 << " ";
                fp_pre << (double)l_pre(0) / lrskid.n_elem * 100 << " " << (double)l_pre(1) / lrskid.n_elem * 100 << " "
                       << (double)l_pre(2) / lrskid.n_elem * 100 << " " << (double)l_pre(3) / lrskid.n_elem * 100 << " ";
                fp_pre << (double)all_pre(0) / N * 100 << " " << (double)all_pre(1) / N * 100 << " " << (double)all_pre(2) / N * 100 << " " << (double)all_pre(3) / N * 100;

                fp_pre << " " << (double)h_AMR(0) / hrskid.n_elem * 100 << " " << (double)h_AMR(1) / hrskid.n_elem * 100 << " "
                       << (double)h_AMR(2) / hrskid.n_elem * 100 << " " << (double)h_AMR(3) / hrskid.n_elem * 100 << " ";
                fp_pre << (double)l_AMR(0) / lrskid.n_elem * 100 << " " << (double)l_AMR(1) / lrskid.n_elem * 100 << " "
                       << (double)l_AMR(2) / lrskid.n_elem * 100 << " " << (double)l_AMR(3) / lrskid.n_elem * 100 << " ";
                fp_pre << (double)all_AMR(0) / N * 100 << " " << (double)all_AMR(1) / N * 100 << " " << (double)all_AMR(2) / N * 100 << " " << (double)all_AMR(3) / N * 100 << " ";

                fp_pre << (int)total_AMR(0) << " " << (int)total_AMR(1) << " " << (int)total_AMR(2) << " " << (int)total_AMR(3) << std::endl;
            }
            //calculate incidece and print out
            if (d % 360 == 0)
            {
                if (fp_ind.is_open())
                {
                    h_ind = sum(diagnosis.rows(hrskid), 0);
                    l_ind = sum(diagnosis.rows(lrskid), 0);
                    all_ind = sum(diagnosis, 0);
                    fp_ind << (int)d << " " << (double)h_ind(0) / hrskid.n_elem * 100 << " " << (double)h_ind(1) / hrskid.n_elem * 100 << " " << (double)h_ind(2) / hrskid.n_elem * 100 << " " << (double)h_ind(3) / hrskid.n_elem * 100 << " ";
                    fp_ind << (double)l_ind(0) / lrskid.n_elem * 100 << " " << (double)l_ind(1) / lrskid.n_elem * 100 << " " << (double)l_ind(2) / lrskid.n_elem * 100 << " " << (double)l_ind(3) / lrskid.n_elem * 100 << " ";
                    fp_ind << (double)all_ind(0) / N * 100 << " " << (double)all_ind(1) / N * 100 << " " << (double)all_ind(2) / N * 100 << " " << (double)all_ind(3) / N * 100;

                    h_ind = sum(diagnosis_AMR.rows(hrskid), 0);
                    l_ind = sum(diagnosis_AMR.rows(lrskid), 0);
                    all_ind = sum(diagnosis_AMR, 0);
                    fp_ind << " " << (double)h_ind(0) / hrskid.n_elem * 100 << " " << (double)h_ind(1) / hrskid.n_elem * 100 << " " << (double)h_ind(2) / hrskid.n_elem * 100 << " " << (double)h_ind(3) / hrskid.n_elem * 100 << " ";
                    fp_ind << (double)l_ind(0) / lrskid.n_elem * 100 << " " << (double)l_ind(1) / lrskid.n_elem * 100 << " " << (double)l_ind(2) / lrskid.n_elem * 100 << " " << (double)l_ind(3) / lrskid.n_elem * 100 << " ";
                    fp_ind << (double)all_ind(0) / N * 100 << " " << (double)all_ind(1) / N * 100 << " " << (double)all_ind(2) / N * 100 << " " << (double)all_ind(3) / N * 100;

                    ind_high = mean(incd_num_test.rows(hrskid), 0) * 100;
                    //ind_high.insert_rows(ind_high.n_rows,mean(incd_num_test.rows(hrskid),0)*100);
                    fp_ind << "   " << (double)ind_high(0) << " " << (double)ind_high(1) << " " << (double)ind_high(2) << " " << (double)ind_high(3) << std::endl;
                }
                diagnosis.zeros();
                diagnosis_AMR.zeros();
                incd_num_test.zeros();
                id_test.ones();
            }
            if (d > par_day + back_years * 360)
            {
                outbreak_d_temp(0) = all_AMR(3);
                outbreak_d_temp(1) = total_AMR(3);
                outbreak_d.insert_cols(outbreak_d.n_cols, outbreak_d_temp);
                if (d - par_day - back_years * 360 == 360 * my_years || all(all_AMR == 0))
                {
                    output(1) = d - par_day - back_years * 360;
                    break;
                }
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "Calculation of indicence  and prevalence: End" << std::endl;

        //    std::cout << "Screening: Begin" << std::endl;
        //}
        //screening
        uvec Scrid = find(scrday == d);
        if (!Scrid.is_empty())
        {
            for (unsigned i = 0; i < Scrid.n_elem; i++)
            {
                uvec ps = find(men[Scrid(i)].status > 0);
                Num_tt(1) = Num_tt(1) + 1; //testing times
                //id_test = 0 not be tested; =1, had a negative test; =2 had a postive test;
                last_testday(Scrid(i)) = d;
                incd_lasttest(Scrid(i), 0) = men[Scrid(i)].status(0) > 0;
                incd_lasttest(Scrid(i), 1) = men[Scrid(i)].status(1) > 0;
                incd_lasttest(Scrid(i), 2) = men[Scrid(i)].status(2) > 0;
                incd_lasttest(Scrid(i), 3) = any(men[Scrid(i)].status > 0);
                if (id_test(Scrid(i)) == 1)
                {
                    incd_num_test.row(Scrid(i)) = incd_num_test.row(Scrid(i)) + incd_lasttest.row(Scrid(i));
                    id_test(Scrid(i)) = 1 * (ps.n_elem == 0) + 2 * (ps.n_elem > 0);
                }
                else
                {
                    id_test(Scrid(i)) = 1 * (ps.n_elem == 0) + 2 * (ps.n_elem > 0);
                    //incd_num_test.row(Scrid(i))=incd_num_test.row(Scrid(i))+incd_lasttest.row(Scrid(i));
                }
                if (ps.n_elem > 0 && unif01(generator) < par_dignrate) //positive and dignosed
                {
                    //whether the casue will be cultured, will be cultured, but might fail
                    if (any(men[Scrid(i)].resistance > 0))
                    {
                        detectAMR_prob = 1 - pow(1 - culture_rate(0), men[Scrid(i)].resistance(0)) * pow(1 - culture_rate(1), men[Scrid(i)].resistance(1)) * pow(1 - culture_rate(2), men[Scrid(i)].resistance(2));
                        detectAMR = unif01(generator) < detectAMR_prob;
                        if (detectAMR)                                   //detected AMR
                        {                                                //if culture is successful
                            if (!detected1st_treat && !detected1st_test) //only execute once
                            {
                                detected1st_test = true;
                                info_1stdetected(0) = men[Scrid(i)].behaviour;
                                info_1stdetected(1) = men[Scrid(i)].acqrate;
                                info_1stdetected(2) = men[Scrid(i)].testday;
                                info_1stdetected(3) = men[Scrid(i)].status(0);
                                info_1stdetected(4) = men[Scrid(i)].status(1);
                                info_1stdetected(5) = men[Scrid(i)].status(2);
                                info_1stdetected(6) = men[Scrid(i)].resistance(0);
                                info_1stdetected(7) = men[Scrid(i)].resistance(1);
                                info_1stdetected(8) = men[Scrid(i)].resistance(2);
                                info_1stdetected(9) = d - par_day - back_years * 360;
                                all_AMR.zeros();
                                //Implementation of recommened screening
                                if (reconsc)
                                {
                                    for (unsigned n = 0; n < N; n++)
                                    {
                                        if (men[n].acqrate > 9)
                                        {
                                            men[n].testday = 90;
                                        }
                                        else
                                        {
                                            if (men[n].testday > 360)
                                            {
                                                men[n].testday = 360;
                                            }
                                        }
                                    }
                                }
                                //total_AMR.zeros();
                                for (unsigned n = 0; n < N; n++)
                                {
                                    AMR_temp.zeros();
                                    cur_AMR(0) = (men[n].resistance(0) > 0);
                                    cur_AMR(1) = (men[n].resistance(1) > 0);
                                    cur_AMR(2) = (men[n].resistance(2) > 0);
                                    cur_AMR(3) = any(men[n].resistance > 0);
                                    all_AMR = all_AMR + cur_AMR;
                                    AMR_temp = cur_AMR > lastAMR.row(n);
                                    diagnosis_AMR.row(n) = diagnosis_AMR.row(n) + AMR_temp;
                                    //total_AMR = total_AMR + AMR_temp;
                                }
                                for (unsigned s = 0; s < 4; s++)
                                {
                                    info_1stdetected(10 + s) = all_AMR(s);
                                    info_1stdetected(14 + s) = total_AMR(s);
                                }

                                if (men[Scrid(i)].RhisPar.n_rows > 1)
                                {
                                    trace_partners.resize(1);
                                    trace_partners(0) = men[Scrid(i)].RhisPar(men[Scrid(i)].RhisPar.n_rows - 1, 1); //current regular parnter
                                    info_1stdetected(18) = any(men[trace_partners(0)].resistance > 0);
                                }
                                else
                                {
                                    info_1stdetected(18) = 0;
                                }
                            }

                            Num_tt(1) = Num_tt(1) + 1;
                            if (unif01(generator) < par_curerate) //cure resistant sites
                            {

                                for (unsigned j = 0; j < 3; j++)
                                {
                                    if (men[Scrid(i)].status(j) > 0)
                                    {
                                        men[Scrid(i)].status(j) = -1;
                                        men[Scrid(i)].periods(j, 3) = d + 7;
                                        men[Scrid(i)].periods(j, 2) = 0;
                                        men[Scrid(i)].periods(j, 1) = 0;
                                        men[Scrid(i)].periods(j, 0) = 0;
                                        men[Scrid(i)].resistance.zeros();
                                    }
                                }
                            }
                            //it's under detectAMR (contact tracing)
                            trace_partners.resize(1);
                            if (men[Scrid(i)].behaviour != 2 && men[Scrid(i)].RhisPar.n_rows > 1 && regular)
                            {
                                trace_partners(0) = men[Scrid(i)].RhisPar(men[Scrid(i)].RhisPar.n_rows - 1, 1); //current regular parnter
                                if (unif01(generator) < par_rr && last_testday(trace_partners(0)) + 14 < d)
                                {
                                    scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //regular partner come to test
                                }
                            }
                            if (men[Scrid(i)].behaviour != 1 && men[Scrid(i)].ChisPar.n_rows > 1 && casual)
                            {
                                for (unsigned k = 0; k < men[Scrid(i)].ChisPar.n_rows - 1; k++)
                                {
                                    if (men[Scrid(i)].ChisPar(men[Scrid(i)].ChisPar.n_rows - 1 - k, 3) > d - 60)
                                    {
                                        trace_partners(0) = men[Scrid(i)].ChisPar(men[Scrid(i)].ChisPar.n_rows - 1 - k, 1); //current casual parnter
                                        if (unif01(generator) < par_rcr && last_testday(trace_partners(0)) + 14 < d)
                                        {
                                            scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //casual partners come to test
                                        }
                                    }
                                }
                            }
                            // group sex
                            if (men[Scrid(i)].behaviour != 1 && men[Scrid(i)].GhisPar.n_rows > 1 && casual && gs_trace)
                            {
                                for (unsigned k = 0; k < men[Scrid(i)].GhisPar.n_rows - 1; k++)
                                {
                                    if (men[Scrid(i)].GhisPar(men[Scrid(i)].GhisPar.n_rows - 1 - k, 3) > d - 60)
                                    {
                                        trace_partners(0) = men[Scrid(i)].GhisPar(men[Scrid(i)].GhisPar.n_rows - 1 - k, 1);
                                        if (unif01(generator) < par_rcr && last_testday(trace_partners(0)) + 14 < d)
                                        {
                                            scrday(trace_partners(0)) = d + ceil(unif01(generator) * 13); //come to test
                                        }
                                    }
                                }
                            }
                            //group sex
                        }
                        else //have resistance but didnot send to curltuere
                        {
                            Num_tt(1) = Num_tt(1) + 1;
                            if (unif01(generator) < par_curerate) //cured non-resistant site
                            {
                                for (unsigned j = 0; j < 3; j++)
                                {
                                    if (men[Scrid(i)].status(j) > 0 && men[Scrid(i)].resistance(j) == 0)
                                    {
                                        men[Scrid(i)].status(j) = -1;
                                        men[Scrid(i)].periods(j, 3) = d + 7;
                                        men[Scrid(i)].periods(j, 2) = 0;
                                        men[Scrid(i)].periods(j, 1) = 0;
                                        men[Scrid(i)].periods(j, 0) = 0;
                                    }
                                }
                            }
                        }
                    }
                    else //no resistance sites
                    {

                        Num_tt(1) = Num_tt(1) + 1;
                        if (unif01(generator) < par_curerate) //cured
                        {
                            for (unsigned j = 0; j < 3; j++)
                            {
                                if (men[Scrid(i)].status(j) > 0)
                                {
                                    men[Scrid(i)].status(j) = -1;
                                    men[Scrid(i)].periods(j, 3) = d + 7;
                                    men[Scrid(i)].periods(j, 2) = 0;
                                    men[Scrid(i)].periods(j, 1) = 0;
                                    men[Scrid(i)].periods(j, 0) = 0;
                                }
                            }
                        }
                    }
                }
                scrday(Scrid(i)) = d + ceil((rnorm(generator) / 10 + 1) * men[Scrid(i)].testday);
            }
        }

        //if (Debuge)
        //{
        //    std::cout << "Screening: End" << std::endl;
        //}
        //introduce backgroud infectious after the build-in period

        if (d == par_day && background)
        {
            if (false)
            {
                std::cout << "Introduce infectious: Begin" << std::endl;
                cout << unif01(generator) << endl;
                cout << randu<double>() << endl;
            }
            uvec siteN = {700, 50, 700};
            vec as_s = {0, 0.9, 0.12};
            unsigned pd;
            int sympt;
            unsigned pd0;
            uvec introInfec = find(randu<vec>(N, 1) < (double)siteN(0) / N);
            for (unsigned j = 0; j < introInfec.n_elem; j++)
            {
                men[introInfec(j)].status(0) = 2;
                laststatus(introInfec(j), 0) = 1;
                laststatus(introInfec(j), 3) = 1;
                pd = gO(generator);
                pd0 = ceil((pd + 3) * unif01(generator));
                men[introInfec(j)].periods(0, 0) = d + 4 - pd0;
                men[introInfec(j)].periods(0, 1) = d + 4 + pd - pd0;
                men[introInfec(j)].periods(0, 2) = d + 4 + pd - pd0;

                men[introInfec(j)].periods(0, 3) = d + 11 + pd - pd0;
            }
            introInfec = find(randu<vec>(N, 1) < (double)siteN(1) / N);
            for (unsigned j = 0; j < introInfec.n_elem; j++)
            {
                sympt = 1 + (unif01(generator) > 0.9);
                men[introInfec(j)].status(1) = sympt;
                laststatus(introInfec(j), 1) = 1;
                laststatus(introInfec(j), 3) = 1;
                pd = gS(generator) * (sympt == 1) + gAU(generator) * (sympt == 2);
                pd0 = ceil((pd + 3) * unif01(generator));
                men[introInfec(j)].periods(1, 0) = d + 4 - pd0;
                men[introInfec(j)].periods(1, 1) = d + 4 + pd - pd0;
                men[introInfec(j)].periods(1, 2) = d + 4 + pd - pd0;

                men[introInfec(j)].periods(1, 3) = d + 11 + pd - pd0;
            }
            introInfec = find(randu<vec>(N, 1) < (double)siteN(2) / N);
            for (unsigned j = 0; j < introInfec.n_elem; j++)
            {
                sympt = 1 + (unif01(generator) > 0.12);
                men[introInfec(j)].status(2) = sympt;
                laststatus(introInfec(j), 2) = 1;
                laststatus(introInfec(j), 3) = 1;
                pd = gS(generator) * (sympt == 1) + gR(generator) * (sympt == 2);
                pd0 = ceil((pd + 3) * unif01(generator));
                men[introInfec(j)].periods(2, 0) = d + 4 - pd0;
                men[introInfec(j)].periods(2, 1) = d + 4 + pd - pd0;
                men[introInfec(j)].periods(2, 2) = d + 4 + pd - pd0;

                men[introInfec(j)].periods(2, 3) = d + 11 + pd - pd0;
            }
            //if (Debuge)
            //{
            //    std::cout << "Introduce infectious: End" << std::endl;
            //}
        }
        //introduce first resistant case,
        if (d == par_day + back_years * 360)
        {
            uvec phi, urt, rec;
            for (unsigned i = 0; i < N; i++)
            {
                if (men[i].status(0) > 0)
                {
                    uvec_push(phi, i);
                }
                if (men[i].status(1) > 0)
                {
                    uvec_push(urt, i);
                }
                if (men[i].status(2) > 0)
                {
                    uvec_push(rec, i);
                }
            }

            unsigned site_1st;
            unsigned pd;
            int sympt;
            unsigned pd0;
            U = unif01(generator);
            site_1st = 0 * (U <= 0.52) + 1 * (U > 0.52 & U <= 0.73) + 2 * (U > 0.73);
            uword infect_1st;
            if (site_1st == 0)
            {
                infect_1st = phi(floor(unif01(generator) * phi.n_elem));
                men[infect_1st].resistance(0) = 1;
            }
            else if (site_1st == 1)
            {

                infect_1st = urt(floor(unif01(generator) * urt.n_elem));
                men[infect_1st].resistance(1) = 1;
            }
            else if (site_1st == 2)
            {
                infect_1st = rec(floor(unif01(generator) * rec.n_elem));
                men[infect_1st].resistance(2) = 1;
            }
            //for other site of selected person
            if (men[infect_1st].status(0) > 0)
            {
                men[infect_1st].resistance(0) = 1;
            }
            if (men[infect_1st].status(1) > 0)
            {
                sympt = 1 + (unif01(generator) > 0.9);
                men[infect_1st].status(1) = sympt;
                pd = gS(generator) * (sympt == 1) + gAU(generator) * (sympt == 2);
                pd0 = ceil((pd + 3) * unif01(generator));
                if (sympt == 1)
                {
                    U = unif01(generator);
                    men[infect_1st].periods(1, 0) = d + ceil(-log(0.65 * U) * 4 * unif01(generator)) - pd0;
                    men[infect_1st].periods(1, 1) = d + ceil(-log(0.65 * U) * 4) - pd0;
                    men[infect_1st].periods(1, 2) = d + ceil(-log(0.65 * U) * 4) + pd - pd0;
                    men[infect_1st].periods(1, 3) = d + ceil(-log(0.65 * U) * 4) + pd + 7 - pd0;
                }
                else if (sympt == 2)
                {
                    men[infect_1st].periods(1, 0) = d + 4 - pd0;
                    men[infect_1st].periods(1, 1) = d + 4 + pd - pd0;
                    men[infect_1st].periods(1, 2) = d + 4 + pd - pd0;
                    men[infect_1st].periods(1, 3) = d + 11 + pd - pd0;
                }
            }
            if (men[infect_1st].status(2) > 0)
            {
                sympt = 1 + (unif01(generator) > 0.12);
                men[infect_1st].status(2) = sympt;
                pd = gS(generator) * (sympt == 1) + gR(generator) * (sympt == 2);
                pd0 = ceil((pd + 3) * unif01(generator));

                if (sympt == 1)
                {
                    U = unif01(generator);
                    men[infect_1st].periods(2, 0) = d + ceil(-log(0.65 * U) * 4 * unif01(generator)) - pd0;
                    men[infect_1st].periods(2, 1) = d + ceil(-log(0.65 * U) * 4) - pd0;
                    men[infect_1st].periods(2, 2) = d + ceil(-log(0.65 * U) * 4) + pd - pd0;
                    men[infect_1st].periods(2, 3) = d + ceil(-log(0.65 * U) * 4) + pd + 7 - pd0;
                }
                else if (sympt == 2)
                {
                    men[infect_1st].periods(2, 0) = d + 4 - pd0;
                    men[infect_1st].periods(2, 1) = d + 4 + pd - pd0;
                    men[infect_1st].periods(2, 2) = d + 4 + pd - pd0;
                    men[infect_1st].periods(2, 3) = d + 11 + pd - pd0;
                }
            }
            //men[infect_1st].status.print();

            trans_tree(0, 1) = infect_1st;
            //have deleted some code
            for (unsigned i = 0; i < N; i++)
            {
                //Implementation of recommened screening

                if (i == infect_1st)
                    continue;
                men[i].status(0) = 0;
                men[i].status(1) = 0;
                men[i].status(2) = 0;
            }

            info_1stinfection(0) = men[infect_1st].behaviour;
            info_1stinfection(1) = men[infect_1st].acqrate;
            info_1stinfection(2) = men[infect_1st].testday;
            info_1stinfection(3) = men[infect_1st].status(0);
            info_1stinfection(4) = men[infect_1st].status(1);
            info_1stinfection(5) = men[infect_1st].status(2);
            info_1stinfection(6) = site_1stAMR;
            //men[infect_1st].RhisPar.print("regular partnerships");
        }
        //Update laststatus
        for (unsigned i = 0; i < N; i++)
        {
            laststatus(i, 0) = (men[i].status(0) > 0);
            laststatus(i, 1) = (men[i].status(1) > 0);
            laststatus(i, 2) = (men[i].status(2) > 0);
            laststatus(i, 3) = any(men[i].status > 0);
            lastAMR(i, 0) = (men[i].resistance(0) > 0);
            lastAMR(i, 1) = (men[i].resistance(1) > 0);
            lastAMR(i, 2) = (men[i].resistance(2) > 0);
            lastAMR(i, 3) = any(men[i].resistance > 0);
        }
    }

    output(0) = runid;
    output.insert_rows(2, info_1stinfection);
    output.insert_rows(output.n_rows, info_1stdetected);
    output.insert_rows(output.n_rows, Num_tt);

    output_infected(0) = runid; //AMR infected
    output_infected.insert_rows(1, all_AMR.t());
    output_infected.insert_rows(output_infected.n_rows, total_AMR.t());
    output_infected.insert_rows(output_infected.n_rows, total_trans_sites.t());

    if (fp_pre.is_open())
    {
        fp_pre.close();
    }
    if (fp_ind.is_open())
    {
        fp_ind.close();
    }
    return;
    delete men;
}

double GPareto_acqrate(int k, int sigma, int theta, long long rseed)
{
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double U = 1 - unif01(generator) * 0.2761495; //(0,0.28)
    double X;
    X = theta + (double)sigma * (pow(U, -k) - 1) / k;
    while (X > 60 || X < 1)
    {
        U = 1 - unif01(generator) * 0.2761495; //(0,0.28)
        X = theta + (double)sigma * (pow(U, -k) - 1) / k;
    }
    return X;
}
Nodeman::Nodeman()
{
    this->name = 0;
}

Nodeman::Nodeman(unsigned int id, long long rseed)
{

    default_random_engine generator(rseed);
    uniform_int_distribution<long long> rseed1(1, 4294967296);
    this->name = id;
    std::uniform_int_distribution<unsigned int> unifage(5400, 23400);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    this->age = unifage(generator);
    double rcb[3] = {0.27, 0.682, 1}; //{0.27,0.4120,0.318};//regular,casual+single,both
    //declaration of acqurate, behaviour

    double U = unif01(generator);
    if (U < rcb[0])
    {
        this->behaviour = 1;
        this->acqrate = 0;

        double anal_condom1[2] = {0.0861, 1}; //{consistely use, sometimes not}
        U = unif01(generator);
        this->anal_condom = 1 - 1 * (U < anal_condom1[0]);
    }
    else if (U < rcb[1] && U >= rcb[0])
    {
        this->behaviour = 2;
        U = unif01(generator);
        if (U < 0.4)
        {
            this->acqrate = unif01(generator) * 0.5;
        }
        else if (U >= 0.4 && U < 1)
        {
            this->acqrate = GPareto_acqrate(15, 7, 1, rseed1(generator));
        }
        double anal_condom2[2] = {0.3712, 1}; //{consistely use, sometimes not}
        U = unif01(generator);
        this->anal_condom = 1 - 1 * (U < anal_condom2[0]);
    }
    else if (U < rcb[2] && U >= rcb[1])
    {
        this->behaviour = 3;
        this->acqrate = GPareto_acqrate(15, 7, 1, rseed1(generator));
        double anal_condom3[2] = {0.3712, 1}; //{consistely use, sometimes not}
        U = unif01(generator);
        this->anal_condom = 1 - 1 * (U < anal_condom3[0]);
    }

    double analdist[3] = {0.11, 0.30, 1}; // The Distribution of anal genders
    U = unif01(generator);
    this->analrole = 1 * (U < analdist[0]) + 2 * (U >= analdist[0] && U < analdist[1]) +
                     3 * (U >= analdist[1] && U <= analdist[2]);

    //test frequency
    if (this->acqrate > 9)
    {
        double testdist_hr[4] = {0.05, 0.49, 0.71, 1.0};
        U = unif01(generator);
        this->testday = (360 + 360 * ceil(unif01(generator))) * (U < 0.05) + 360 * (U >= 0.05 && U < 0.49) + 180 * (U >= 0.49 && U < 0.71) + 120 * (U >= 0.71 && U < 0.89) + 90 * (U >= 0.89 && U < 1);
    }
    else if (this->acqrate > 0 && this->acqrate < 1)
    {
        this->testday = (360 + 360 * ceil(unif01(generator)));
    }
    else
    {
        this->testday = (360 + 360 * ceil(unif01(generator))) * (U < 0.10) + 360 * (U >= 0.10 && U < 1);
    }
}

void Nodeman::Initializer(unsigned int id, long long rseed)
{
    default_random_engine generator(rseed);
    uniform_int_distribution<long long> rseed1(1, 4294967296);
    name = id;
    std::uniform_int_distribution<unsigned int> unifage(5400, 23400);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    age = unifage(generator);

    double rcb[3] = {0.28, 0.61, 1}; //{0.28,0.40,0.32};//regular,casual+single,both
    //declaration of acqurate, behaviour
    double U = unif01(generator);
    if (U < rcb[0])
    {
        behaviour = 1;
        acqrate = 0;

        double anal_condom1[2] = {0.0861, 1}; //{consistely use 1, sometimes not 0}
        U = unif01(generator);
        anal_condom = (U < anal_condom1[0]);
    }
    else if (U < rcb[1] && U >= rcb[0])
    {
        behaviour = 2;
        U = unif01(generator);
        if (U < 0.34)
        {
            acqrate = unif01(generator) * 0.5;
        }
        else if (U >= 0.34 && U < 1)
        {
            acqrate = GPareto_acqrate(15, 7, 1, rseed1(generator));
        }
        double anal_condom2[2] = {0.3712, 1}; //{consistely use 1, sometimes not 0}
        U = unif01(generator);
        anal_condom = 1 * (U < anal_condom2[0]);
    }
    else if (U < rcb[2] && U >= rcb[1])
    {
        behaviour = 3;
        acqrate = GPareto_acqrate(15, 7, 1, rseed1(generator));
        double anal_condom3[2] = {0.3712, 1}; //{consistely use 1, sometimes not 1}
        U = unif01(generator);
        anal_condom = 1 * (U < anal_condom3[0]);
    }

    double analdist[3] = {0.11, 0.30, 1}; // The Distribution of anal genders
    U = unif01(generator);
    analrole = 1 * (U < analdist[0]) + 2 * (U >= analdist[0] && U < analdist[1]) +
               3 * (U >= analdist[1] && U <= analdist[2]);

    //test frequency
    if (this->acqrate > 9)
    {
        double testdist_hr[4] = {0.05, 0.49, 0.71, 1.0};
        U = unif01(generator);
        //testday = (360)*(U<0.05)
        //+360*(U>=0.05 && U<0.49) +180*(U>=0.49 && U<0.71)
        //+120*(U>=0.71 && U<0.89) + 90*(U>=0.89 && U<1);
        testday = (480 + 360 * ceil(unif01(generator))) * (U < 0.05) + 360 * (U >= 0.05 && U < 0.49) + 180 * (U >= 0.49 && U < 0.71) + 120 * (U >= 0.71 && U < 0.89) + 90 * (U >= 0.89 && U < 1);
    }
    else if (this->acqrate > 0 && this->acqrate < 1)
    {
        testday = (480 + 360 * ceil(unif01(generator)));
    }
    else
    {
        testday = (360 + 360 * ceil(unif01(generator))) * (U < 0.25) //was 0.10
                  + 360 * (U >= 0.25 && U < 1);
    }
}

unsigned int getstopday(unsigned int t1, unsigned int t2, unsigned int t3)
{
    unsigned int minvalue = t1;
    if (minvalue > t2)
    {
        minvalue = t2;
    }
    if (minvalue > t3)
    {
        minvalue = t3;
    }
    return minvalue;
}
void uvec_push(arma::uvec &v, unsigned int value)
{
    arma::uvec av(1);
    av.at(0) = value;
    v.insert_rows(v.n_rows, av.row(0));
}
void umat_push(arma::umat &v, unsigned int value1, unsigned int value2)
{
    arma::uvec av(2);
    av.at(0) = value1;
    av.at(1) = value2;
    v.insert_rows(v.n_rows, av.t());
}

void rirj_partnership(Nodeman &mri, Nodeman &mrj, mat &eventsij, unsigned int &stopDayR, uword &ratri, uword &ratrj, unsigned int today,
                      long long rseed)
{
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    unsigned int stopd;
    urowvec temPartner(4);
    mat events(5, 2);
    events.zeros();
    //typeM = 1;
    stopd = ceil(-log(unif01(generator)) * 1440) + 1;

    stopDayR = today + stopd; //StopDayR(ri, rj)
    ratri = today + stopd;    //Rat(ri) = stopd;
    ratrj = today + stopd;    //Rat(rj) = stopd;
    //temPartner = {mri.name, mrj.name, today, stopDayR};
    temPartner(0) = mri.name;
    temPartner(1) = mrj.name;
    temPartner(2) = today;
    temPartner(3) = stopDayR;
    mri.RhisPar.insert_rows(mri.RhisPar.n_rows, temPartner);
    //temPartner = {mrj.name, mri.name, today, stopDayR};
    temPartner(0) = mrj.name;
    temPartner(1) = mri.name;
    temPartner(2) = today;
    temPartner(3) = stopDayR;
    mrj.RhisPar.insert_rows(mrj.RhisPar.n_rows, temPartner);
    if (unif01(generator) < 0.83)
    {
        events(0, 0) = (double)(2.4 + 1.2 * unif01(generator)) / 7;
    };
    if (unif01(generator) < 0.825)
    {
        events(1, 0) = (double)(1.6 + 0.6 * unif01(generator)) / 7;
    };
    if (unif01(generator) < 0.60)
    {
        events(2, 0) = (double)(1.2 + 0.6 * unif01(generator)) / 7;
    };
    double prob_anal;
    if (mri.analrole == 1 && mrj.analrole == 1)
    {
        prob_anal = 0;
    }
    else if (mri.analrole == 2 && mrj.analrole == 2)
    {
        prob_anal = 0;
    }
    else
    {
        prob_anal = 0.50;
    };
    if (unif01(generator) < prob_anal)
    {
        events(3, 0) = (double)(1.6 + 0.8 * unif01(generator)) / 7;
    };
    if (unif01(generator) < 0.603)
    {
        events(4, 0) = (double)(1.2 + 0.6 * unif01(generator)) / 7;
    };
    eventsij = events;
}
void cicj_partnership(Nodeman &mci, Nodeman &mcj, mat &eventsij, unsigned int &stopDayC, uword &catci, uword &catcj, unsigned int today,
                      long long rseed)
{
    //typeM = 2; //TypeM(ci,cj)=2;
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    catci = ceil(-log(unif01(generator)) * (double)180 / mci.acqrate) + today;
    catcj = ceil(-log(unif01(generator)) * (double)180 / mcj.acqrate) + today;
    unsigned int sj_overall = ceil(unif01(generator) * 14); //ceil(-log(unif01(generator)) * 14) ;
    unsigned int stopdaycsj = getstopday(catci, catcj, sj_overall + today) - today;
    stopDayC = stopdaycsj + today; //StopDayC(ci,cj)
    urowvec temPartner(4);
    mat events(5, 2);
    events.zeros();
    //temPartner = {mci.name, mcj.name, today, stopDayC};
    temPartner(0) = mci.name;
    temPartner(1) = mcj.name;
    temPartner(2) = today;
    temPartner(3) = stopDayC;
    mci.ChisPar.insert_rows(mci.ChisPar.n_rows, temPartner);
    //temPartner = {mcj.name, mci.name, today, stopDayC};
    temPartner(0) = mcj.name;
    temPartner(1) = mci.name;
    temPartner(2) = today;
    temPartner(3) = stopDayC;
    mcj.ChisPar.insert_rows(mcj.ChisPar.n_rows, temPartner);
    double pday;
    double prd;

    if (unif01(generator) < 0.83)
    {
        pday = (double)(2.4 + 1.2 * unif01(generator)) / 7;
        prd = (double)(pday * stopdaycsj - 1) / (stopdaycsj - 1);
        events(0, 0) = (prd > 0) * prd + (prd < 0) * 0;
        events(0, 1) = 1;
    };
    if (unif01(generator) < 0.825)
    {
        pday = (double)(1.6 + 0.6 * unif01(generator)) / 7;
        prd = (double)(pday * stopdaycsj - 1) / (stopdaycsj - 1);
        events(1, 0) = (prd > 0) * prd + (prd < 0) * 0;
        events(1, 1) = 1;
    };
    if (unif01(generator) < 0.60)
    {
        pday = (double)(1.2 + 0.6 * unif01(generator)) / 7;
        prd = (double)(pday * stopdaycsj - 1) / (stopdaycsj - 1);
        events(2, 0) = (prd > 0) * prd + (prd < 0) * 0;
        events(2, 1) = 1;
    };
    double prob_anal;
    if (mci.analrole == 1 && mcj.analrole == 1)
    {
        prob_anal = 0;
    }
    else if (mci.analrole == 2 && mcj.analrole == 2)
    {
        prob_anal = 0;
    }
    else
    {
        prob_anal = 0.50;
    };
    if (unif01(generator) < prob_anal)
    {
        pday = (double)(1.6 + 0.8 * unif01(generator)) / 7;
        prd = (double)(pday * stopdaycsj - 1) / (stopdaycsj - 1);
        events(3, 0) = (prd > 0) * prd + (prd < 0) * 0;
        events(3, 1) = 1;
    };
    if (unif01(generator) < 0.603)
    {
        pday = (double)(1.2 + 0.6 * unif01(generator)) / 7;
        prd = (double)(pday * stopdaycsj - 1) / (stopdaycsj - 1);
        events(4, 0) = (prd > 0) * prd + (prd < 0) * 0;
        events(4, 1) = 1;
    };
    eventsij = events; //EvM(ci,cj)
}

void gigj_partnership(Nodeman &mgi, Nodeman &mgj, mat &eventsij, rowvec &pgi_ac, rowvec &pgj_ac, unsigned int today,
                      long long rseed)
{
    //typeM = 3; //TypeM(gi,gj)=3;
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    urowvec temPartner;
    mat events(5, 2);
    events.zeros();
    temPartner = {mgi.name, mgj.name, today + 1, today + 1};
    mgi.GhisPar.insert_rows(mgi.GhisPar.n_rows, temPartner);
    temPartner = {mgj.name, mgi.name, today + 1, today + 1};
    mgj.GhisPar.insert_rows(mgj.GhisPar.n_rows, temPartner);
    double pday;
    double prd;
    if (unif01(generator) < 0.83)
    {
        events(0, 1) = 1;
    };
    if (unif01(generator) < 0.825)
    {
        events(1, 1) = 1;
    };
    if (unif01(generator) < 0.60)
    {
        events(2, 1) = 1;
    }
    double prob_anal;
    if (mgi.analrole == 1 && mgj.analrole == 1)
    {
        prob_anal = 0;
    }
    else if (mgi.analrole == 2 && mgj.analrole == 2)
    {
        prob_anal = 0;
    }
    else
    {
        prob_anal = 0.50;
    }
    if (unif01(generator) < prob_anal)
    {
        events(3, 1) = 1;
    }
    if (unif01(generator) < 0.603)
    {
        events(4, 1) = 1;
    }
    eventsij = events; //EvM(gi,gj)

    unsigned NumPar6m = mgi.ChisPar.n_rows + mgi.GhisPar.n_rows - 2;
    if (NumPar6m > mgi.acqrate)
    {
        pgi_ac(0) = 0;
    }
    else
    {
        unsigned dayC = 0, dayG = 0;
        if (mgi.ChisPar.n_rows > 1)
        {
            dayC = mgi.ChisPar(1, 3);
        }
        if (mgi.GhisPar.n_rows > 1)
        {
            dayG = mgi.GhisPar(1, 3);
        }
        if (dayC > 0 && dayG == 0)
        {
            pgi_ac(0) = (double)(mgi.acqrate - NumPar6m) / (180 - today + dayC);
        }
        else if (dayC == 0 && dayG > 0)
        {
            pgi_ac(0) = (double)(mgi.acqrate - NumPar6m) / (180 - today + dayG);
        }
        else if (dayC > 0 && dayG > 0)
        {
            pgi_ac(0) = (double)(mgi.acqrate - NumPar6m) / (180 - today + min({dayC, dayG}));
        }
        else
        {
            pgi_ac(0) = (double)1 / (180 / mgi.acqrate - today + pgi_ac(1));
            pgi_ac(0) = 1 * (pgi_ac(0) < 0) + pgi_ac(0) * (pgi_ac(0) > 0);
            //why here is 1?
        }
    }
    NumPar6m = mgj.ChisPar.n_rows + mgj.GhisPar.n_rows - 2;
    if (NumPar6m > mgj.acqrate)
    {
        pgj_ac(0) = 0;
    }
    else
    {
        unsigned dayC = 0, dayG = 0;
        if (mgj.ChisPar.n_rows > 1)
        {
            dayC = mgj.ChisPar(1, 3);
        }
        if (mgj.GhisPar.n_rows > 1)
        {
            dayG = mgj.GhisPar(1, 3);
        }
        if (dayC > 0 && dayG == 0)
        {
            pgj_ac(0) = (double)(mgj.acqrate - NumPar6m) / (180 - today + dayC);
        }
        else if (dayC == 0 && dayG > 0)
        {
            pgj_ac(0) = (double)(mgj.acqrate - NumPar6m) / (180 - today + dayG);
        }
        else if (dayC > 0 && dayG > 0)
        {
            pgj_ac(0) = (double)(mgj.acqrate - NumPar6m) / (180 - today + min({dayC, dayG}));
        }
        else
        {
            pgj_ac(0) = (double)1 / ((double)180 / mgi.acqrate - today + pgj_ac(1));
            pgj_ac(0) = 1 * (pgj_ac(0) < 0) + pgj_ac(0) * (pgj_ac(0) > 0);
            //why here is 1?
        }
    }
}

mat sch_sex(mat events, long long rseed)
{
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    events.col(1).zeros();
    if (events(0, 0) > unif01(generator))
    {
        events(0, 1) = 1;
    }
    if (events(1, 0) > unif01(generator))
    {
        events(1, 1) = 1;
    }
    if (events(2, 0) > unif01(generator))
    {
        events(2, 1) = 1;
    }
    if (events(3, 0) > unif01(generator))
    {
        events(3, 1) = 1;
    }
    if (events(4, 0) > unif01(generator))
    {
        events(4, 1) = 1;
    }

    return events;
}
//
void transAB(Nodeman &A, Nodeman &B, mat events, unsigned typem, unsigned int today, vec p_act, urowvec &trans_sites,
             long long rseed)
{
    default_random_engine generator(rseed);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    std::uniform_real_distribution<double> gS(1.0, 3.0);
    //std::gamma_distribution<double> gS(3.01, 0.8615); //truncate late.
    std::gamma_distribution<double> gAU(206.4, 0.4118);
    std::gamma_distribution<double> gO(201.6, 0.4167);
    std::gamma_distribution<double> gR(3702.857, 0.0972);

    unsigned int prot_anal;
    unsigned int pd;
    double U;
    if (A.anal_condom == 0 && B.anal_condom == 0)
    {
        prot_anal = unif01(generator) < 0.5;
    }
    else
    {
        prot_anal = 1;
    }

    if (events(0, 1) == 1)
    {
        if (A.status(0) > 0 && B.status(0) == 0)
        { //from A.p to B.p
            if (today > A.periods(0, 0) && today < A.periods(0, 1))
            {
                if (unif01(generator) < p_act(0))
                {
                    B.status(0) = 2;
                    pd = round(gO(generator));
                    B.periods(0, 0) = today + 4;
                    B.periods(0, 1) = today + 4 + pd;
                    B.periods(0, 2) = today + 4 + pd;
                    B.periods(0, 3) = today + 11 + pd;
                    B.resistance(0) = A.resistance(0);
                    trans_sites(0) = 1;
                }
            }
        }
        else if ((A.status(0) == 0 && B.status(0) > 0))
        { //from B.p to A.p
            if (today > B.periods(0, 0) && today < B.periods(0, 1))
            {
                if (unif01(generator) < p_act(0))
                {
                    A.status(0) = 2;
                    pd = round(gO(generator));
                    A.periods(0, 0) = today + 4;
                    A.periods(0, 1) = today + 4 + pd;
                    A.periods(0, 2) = today + 4 + pd;
                    A.periods(0, 3) = today + 11 + pd;
                    A.resistance(0) = B.resistance(0);
                    trans_sites(0) = 1;
                }
            }
        }
    }
    //cout<<"Event11 : End" <<endl;
    if (events(1, 1) == 1)
    {
        if (unif01(generator) < 1)
        { //A is receptive and B is insertive
            if (A.status(0) > 0 && B.status(1) == 0)
            { //A.p to B.u
                if (today > A.periods(0, 0) && today < A.periods(0, 1))
                {
                    if (unif01(generator) < p_act(1))
                    {
                        if (unif01(generator) < 0.9)
                        {
                            B.status(1) = 1;
                            pd = round(gS(generator));
                            U = unif01(generator);
                            B.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            B.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                            B.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            B.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            B.resistance(1) = A.resistance(0);
                        }
                        else
                        {
                            B.status(1) = 2;
                            pd = round(gAU(generator));
                            B.periods(1, 0) = today + 4;
                            B.periods(1, 1) = today + 4 + pd;
                            B.periods(1, 2) = today + 4 + pd;
                            B.periods(1, 3) = today + 11 + pd;
                            B.resistance(1) = A.resistance(0);
                        }
                        trans_sites(1) = 1;
                    }
                }
            }
            else if (A.status(0) == 0 && B.status(1) > 0)
            { // B.u to A.p
                if (today > B.periods(1, 0) && today < B.periods(1, 1))
                {
                    if (unif01(generator) < p_act(2))
                    {
                        A.status(0) = 2;
                        pd = round(gO(generator));
                        A.periods(0, 0) = today + 4;
                        A.periods(0, 1) = today + 4 + pd;
                        A.periods(0, 2) = today + 4 + pd;
                        A.periods(0, 3) = today + 11 + pd;
                        A.resistance(0) = B.resistance(1);
                        trans_sites(2) = 1;
                    }
                }
            }
        }
        if (unif01(generator) < 1) //A is insertive and B is receptive
        {
            if (A.status(1) > 0 && B.status(0) == 0) //from A.u to B.p
            {
                if (today > A.periods(1, 0) && today < A.periods(1, 1))
                {
                    if (unif01(generator) < p_act(2)) //p
                    {                                 // A.u to B.p
                        B.status(0) = 2;
                        pd = round(gO(generator));
                        B.periods(0, 0) = today + 4;
                        B.periods(0, 1) = today + 4 + pd;
                        B.periods(0, 2) = today + 4 + pd;
                        B.periods(0, 3) = today + 11 + pd;
                        B.resistance(0) = A.resistance(1);
                        trans_sites(2) = 1;
                    }
                }
            }
            else if (A.status(1) == 0 && B.status(0) > 1)
            { //from B.p to A.u
                if (today > B.periods(0, 0) && today < B.periods(0, 1))
                {
                    if (unif01(generator) < p_act(1))
                    {
                        if (unif01(generator) < 0.9)
                        {
                            A.status(1) = 1;
                            pd = round(gS(generator));
                            U = unif01(generator);

                            A.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            A.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                            A.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            A.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            A.resistance(1) = B.resistance(0);
                        }
                        else
                        {
                            A.status(1) = 2;
                            pd = round(gAU(generator));
                            A.periods(1, 0) = today + 4;
                            A.periods(1, 1) = today + 4 + pd;
                            A.periods(1, 2) = today + 4 + pd;
                            A.periods(1, 3) = today + 11 + pd;
                            A.resistance(1) = B.resistance(0);
                        }
                        trans_sites(1) = 1;
                    }
                }
            }
        }
    }
    //cout<<"Event12 21: End" <<endl;

    if (events(2, 1) == 1)
    {
        if (unif01(generator) < 1)
        {
            if (A.status(0) > 0 && B.status(2) == 0)
            { //13 from A.p to partner B.r
                if (today > A.periods(0, 0) && today < A.periods(0, 1))
                {
                    if (unif01(generator) < p_act(3))
                    {
                        if (unif01(generator) < 0.12)
                        {
                            B.status(2) = 1;           //periods set and resistance set
                            pd = round(gS(generator)); //ceil(2.6+randn*2.24);
                            U = unif01(generator);
                            B.periods(2, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            B.periods(2, 1) = today + ceil(-log(0.65 * U) * 4);
                            B.periods(2, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            B.periods(2, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            B.resistance(2) = A.resistance(0);
                        }
                        else
                        {
                            B.status(2) = 2;
                            pd = round(gR(generator)); // ceil(336+rand*28);
                            B.periods(2, 0) = today + 4;
                            B.periods(2, 1) = today + pd + 4;
                            B.periods(2, 2) = today + pd + 4;
                            B.periods(2, 3) = today + 11 + pd;
                            B.resistance(2) = A.resistance(0);
                        }
                        trans_sites(3) = 1;
                    }
                }
            }
            else if (A.status(0) == 0 && B.status(2) > 0) ////31 from B.r to partner A.p
            {
                if (today > B.periods(2, 0) && today < B.periods(2, 1))
                {
                    if (unif01(generator) < p_act(4))
                    {
                        A.status(0) = 2;
                        pd = round(gO(generator)); //ceil(70+rand*63);//ceil(84+randn*35);
                        A.periods(0, 0) = today + 4;
                        A.periods(0, 1) = today + pd + 4;
                        A.periods(0, 2) = today + pd + 4;
                        A.periods(0, 3) = today + pd + 11;
                        A.resistance(0) = B.resistance(2);
                        trans_sites(4) = 1;
                    }
                }
            }
        }
        if (unif01(generator) < 1)
        {
            if (A.status(2) > 0 && B.status(0) == 0)
            { //31 from A.r to partner B.p
                if (today > A.periods(2, 0) && today < A.periods(2, 1))
                {
                    if (unif01(generator) < p_act(4))
                    {
                        B.status(0) = 2;
                        pd = round(gO(generator)); //ceil(70+rand*63);//ceil(84+randn*35);
                        B.periods(0, 0) = today + 4;
                        B.periods(0, 1) = today + pd + 4;
                        B.periods(0, 2) = today + pd + 4;

                        B.periods(0, 3) = today + pd + 11;
                        B.resistance(0) = A.resistance(2);
                        trans_sites(4) = 1;
                    }
                }
            }
            else if (A.status(2) == 0 && B.status(0) > 0) //13 from B.p to partner A.r
            {
                if (today > B.periods(0, 0) && today < B.periods(0, 1))
                {
                    if (unif01(generator) < p_act(3))
                    {
                        if (unif01(generator) < 0.12) //
                        {
                            A.status(2) = 1;           //periods set and resistance set
                            pd = round(gS(generator)); //ceil(2.6+randn*2.24);
                            U = unif01(generator);

                            A.periods(2, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            A.periods(2, 1) = today + ceil(-log(0.65 * U) * 4);
                            A.periods(2, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            A.periods(2, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            A.resistance(2) = B.resistance(0);
                        }
                        else
                        {
                            A.status(2) = 2;
                            pd = round(gR(generator)); //ceil(336+rand*28);
                            A.periods(2, 0) = today + 4;
                            A.periods(2, 1) = today + pd + 4;
                            A.periods(2, 2) = today + pd + 4;
                            A.periods(2, 3) = today + pd + 11;
                            A.resistance(2) = B.resistance(0);
                        }
                        trans_sites(3) = 1;
                    }
                }
            }
        }
    }
    //cout<<"Event13 31: End" <<endl;

    if (events(3, 1) == 1) //23 A.u to B.r, A.r to B.u, B.r to A.u, B,u to A.r
    {
        unsigned int anal_pref = 0;
        if (A.analrole == 3 && B.analrole == 3)
        {
            unsigned int U = unif01(generator);
            //if both individuals are versatile, we assume they can do anal intercourse in both ways.
            anal_pref = 1; //1 * (U < 0.5) + 2 * (U > 0.5);
        }

        if ((A.analrole == 1 && (B.analrole == 2 || B.analrole == 3)) || anal_pref == 1)
        { //A is receptive, B is insertive
            //From A r to B u or from B u to A r
            if (A.status(2) > 0 && B.status(1) == 0) //23 A.r to B.u
            {
                if (today > A.periods(2, 0) && today < A.periods(2, 1))
                {
                    if (unif01(generator) < p_act(6) * pow(0.10, prot_anal))
                    {
                        if (unif01(generator) < 0.90) // u
                        {
                            B.status(1) = 1;           //periods set and resistance set
                            pd = round(gS(generator)); //round(6.3+rand*1.4);
                            U = unif01(generator);
                            B.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            B.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                            B.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            B.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            B.resistance(1) = A.resistance(2);
                        }
                        else
                        {
                            B.status(1) = 2;
                            pd = round(gAU(generator)); //round(70+rand*70);
                            B.periods(1, 0) = today + 4;
                            B.periods(1, 1) = today + pd + 4;
                            B.periods(1, 2) = today + pd + 4;

                            B.periods(1, 3) = today + pd + 11;
                            B.resistance(1) = A.resistance(2);
                        }
                        trans_sites(6) = 1;
                    }
                }
            }
            else if (A.status(2) == 0 && B.status(1) > 0) //32 from B.u to A.r
            {
                if (today > B.periods(1, 0) && today < B.periods(1, 1))
                {
                    if (unif01(generator) < p_act(5) * pow(0.10, prot_anal))
                    {
                        if (unif01(generator) < 0.12) //
                        {
                            A.status(2) = 1; //periods set and resistance set
                            pd = round(gS(generator));
                            U = unif01(generator);
                            A.periods(2, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            A.periods(2, 1) = today + ceil(-log(0.65 * U) * 4);
                            A.periods(2, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                            A.periods(2, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            A.resistance(2) = B.resistance(1);
                        }
                        else
                        {
                            A.status(2) = 2;
                            pd = round(gR(generator));
                            A.periods(2, 0) = today + 4;
                            A.periods(2, 1) = today + 4 + pd;

                            A.periods(2, 2) = today + 4 + pd;

                            A.periods(2, 3) = today + 11 + pd;
                            A.resistance(2) = B.resistance(1);
                        }
                        trans_sites(5) = 1;
                    }
                }
            }
        }
        if ((A.analrole == 2 && (B.analrole == 1 || B.analrole == 3)) || anal_pref == 1)
        // A is insertive, B is receptive
        // From A u to B r, or B r to A u
        {
            if (A.status(1) > 0 && B.status(2) == 0) //23 from A.u to partner  B.r
            {
                if (today > A.periods(1, 0) && today < A.periods(1, 1))
                {
                    if (unif01(generator) < p_act(5) * pow(0.10, prot_anal))
                    {
                        if (unif01(generator) < 0.12) //
                        {
                            B.status(2) = 1;           //periods set and resistance set
                            pd = round(gS(generator)); //ceil(2.6+randn*2.24);
                            U = unif01(generator);
                            B.periods(2, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            B.periods(2, 1) = today + ceil(-log(0.65 * U) * 4);
                            B.periods(2, 2) = today + ceil(-log(0.65 * U) * 4) + pd;

                            B.periods(2, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 11;
                            B.resistance(2) = A.resistance(1);
                        }
                        else
                        {
                            B.status(2) = 2;
                            pd = round(gR(generator)); //ceil(336+rand*28);
                            B.periods(2, 0) = today;
                            B.periods(2, 1) = today + 4;
                            B.periods(2, 2) = today + 4;

                            B.periods(2, 3) = today + pd + 11;
                            B.resistance(2) = A.resistance(1);
                        }
                        trans_sites(5) = 1;
                    }
                }
            }
            else if (A.status(1) == 0 && B.status(2) > 0) //32 from B.r to A.u
            {
                if (today > B.periods(2, 0) && today < B.periods(2, 1))
                {
                    if (unif01(generator) < p_act(6) * pow(0.10, prot_anal))
                    {
                        if (unif01(generator) < 0.90) //u
                        {
                            A.status(1) = 1;           //periods set and resistance set
                            pd = round(gS(generator)); //ceil(6.3+rand*1.4);
                            U = unif01(generator);
                            A.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                            A.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                            A.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;

                            A.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                            A.resistance(1) = B.resistance(2);
                        }
                        else
                        {
                            A.status(1) = 2;
                            pd = round(gAU(generator)); //ceil(70+rand*70); //wrong
                            A.periods(1, 0) = today + 4;
                            A.periods(1, 1) = today + pd + 4;
                            A.periods(1, 2) = today + pd + 4;
                            A.periods(1, 3) = today + pd + 11;
                            B.resistance(1) = A.resistance(2);
                        }
                        trans_sites(6) = 1;
                    }
                }
            }
        }
    }
    //cout<<"Event23 32 : End" <<endl;

    if (events(4, 1) == 1) //22
    {
        if (A.status(1) > 0 && B.status(1) == 0) //from A.u to B.u
        {
            if (today > A.periods(1, 0) && today < A.periods(1, 1))
            {
                if (unif01(generator) < p_act(7))
                {
                    if (unif01(generator) < 0.90) // u
                    {
                        B.status(1) = 1;           //periods set and resistance set
                        pd = round(gS(generator)); //round(6.3+rand*1.4);
                        U = unif01(generator);

                        B.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                        B.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                        B.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;

                        B.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                        B.resistance(1) = A.resistance(1);
                    }
                    else
                    {
                        B.status(1) = 2;
                        pd = round(gAU(generator)); //round(70+rand*70);
                        B.periods(1, 0) = today + 4;
                        B.periods(1, 1) = today + pd + 4;
                        B.periods(1, 2) = today + pd + 4;

                        B.periods(1, 3) = today + pd + 11;
                        B.resistance(1) = A.resistance(1);
                    }
                    trans_sites(7) = 1;
                }
            }
        }
        else if (A.status(1) == 0 && B.status(1) > 0) //from B.u to A.u
        {
            if (today > B.periods(1, 0) && today < B.periods(1, 1))
            {
                if (unif01(generator) < p_act(7))
                {

                    if (unif01(generator) < 0.90) //u
                    {
                        A.status(1) = 1;           //periods set and resistance set
                        pd = round(gS(generator)); //ceil(6.3+rand*1.4);
                        U = unif01(generator);

                        A.periods(1, 0) = today + ceil(-log(0.65 * U) * 4 * unif01(generator));
                        A.periods(1, 1) = today + ceil(-log(0.65 * U) * 4);
                        A.periods(1, 2) = today + ceil(-log(0.65 * U) * 4) + pd;
                        A.periods(1, 3) = today + ceil(-log(0.65 * U) * 4) + pd + 7;
                        A.resistance(1) = B.resistance(1);
                    }
                    else
                    {
                        A.status(1) = 2;
                        pd = round(gAU(generator)); //ceil(70+rand*70); //wrong
                        A.periods(1, 0) = today + 4;
                        A.periods(1, 1) = today + pd + 4;
                        A.periods(1, 2) = today + pd + 4;

                        A.periods(1, 3) = today + pd + 11;
                        A.resistance(1) = B.resistance(1);
                    }
                    trans_sites(7) = 1;
                }
            }
        }
    }
    //cout<<"Event11 : End" <<endl;
}
