#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_map>
#include <unordered_set>




using namespace std;




using rng_type  = std::mt19937_64;
using rng_int_t = typename rng_type::result_type;
auto constexpr rng_wsize = rng_type::word_size;


auto file_exist( string const& file_name )
    -> bool
{
    auto ifile = ifstream( file_name );
    return ifile.good();
}

int main()
{
       /**************************************************
        ** priors *****************************************
        **************************************************/
	/// @info : importance of the word in a speech.
	auto L_repr_min = 1;
	auto L_repr_max = 5;

        auto n_Rg     = size_t{1};


        // random number generator. ///////////////////////
        auto rng = rng_type{ random_device{}() };

        //
        uniform_real_distribution<> rng_n_Rl(0,2);
	auto n_Rl = size_t(n_Rg * pow(10, rng_n_Rl(rng))); 

	// population sizes. //////////////////////////////
	uniform_real_distribution<> rng_N0_ind(1,3);
	auto N0 = size_t(pow(10, rng_N0_ind(rng)));
	uniform_real_distribution<> rng_N1_ST_ind(1,4.3);
	auto N1_ST = size_t(pow(10, rng_N1_ST_ind(rng)));
	uniform_real_distribution<> rng_N2_ST_ind(2,5.5);
	auto N2_ST = size_t(pow(10, rng_N2_ST_ind(rng)));
	uniform_real_distribution<> rng_N1_MA_ind(1,3.3);
	auto N1_MA = size_t(pow(10, rng_N1_MA_ind(rng)));
	uniform_real_distribution<> rng_N2_MA_ind(2,4);
	auto N2_MA = size_t(pow(10, rng_N2_MA_ind(rng)));
	uniform_real_distribution<> rng_Nb_MA_ind(1,3);
	auto Nb_MA = size_t(pow(10, rng_Nb_MA_ind(rng)));
	auto N1_ST_MA = N1_ST + N1_MA;
	auto N2_ST_MA = N2_ST + N2_MA;

	// event times. ///////////////////////////////////
	uniform_int_distribution<size_t> rng_t0(815,815);
	auto t0 = rng_t0(rng); 
	uniform_int_distribution<size_t> rng_t_Nexp(15,215);
	auto t_Nexp = rng_t_Nexp(rng);
	uniform_int_distribution<size_t> rng_t_MA(241,815);
	auto t_MA_e = rng_t_MA(rng);
	auto t_MA_b = rng_t_MA(rng);
	while(t_MA_e > t_MA_b)
	{
		t_MA_e = rng_t_MA(rng);
		t_MA_b = rng_t_MA(rng);
	}


	// listening ability. -----------------------------
	uniform_real_distribution<> rng_L_alpha(0.01,1);
	auto L_alpha = rng_L_alpha(rng) ; // weight of the speaker.

	// mutation rates. --------------------------------
	uniform_real_distribution<> rng_mu_l(2,6);
	auto L_mut = double(pow(10, -rng_mu_l(rng)));
	
	uniform_real_distribution<> rng_mu_g(7,9);
	auto G_mut = double( pow(10, -rng_mu_g(rng)) );




auto param_name = "param.txt";
		auto file_existed = file_exist( param_name );
		if( auto param_file = ofstream( param_name, ios::app ); param_file.fail() )
            throw runtime_error( "cannot open parameter file." );
        else
        {
            if( !file_existed )
                param_file << "n_Rg\tn_Rl\tN0\tN1_ST\tN2_ST\tN1_MA\tN2_MA\tNb_MA\tN1_ST_MA\tN2_ST_MA\tt_Nexp\tt_MA_e\tt_MA_b\tt0\tL_alpha\tL_mut\tG_mut\tL_repr_min\tL_repr_max";

            param_file << "\n" << n_Rg << "\t" << n_Rl << "\t" << N0 << "\t" << N1_ST << "\t" << N2_ST << "\t" << N1_MA << "\t" << N2_MA << "\t" <<  Nb_MA  << "\t" << N1_ST_MA  << "\t" << N2_ST_MA  << "\t" << t_Nexp << "\t" << t_MA_e << "\t" << t_MA_b << "\t" << t0  << "\t" << L_alpha << "\t" << L_mut << "\t" << G_mut << "\t" << L_repr_min << "\t" << L_repr_max;
	param_file.close();
}
	return 0;
}

