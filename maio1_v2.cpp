#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <map>




#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/beta.hpp>




using namespace std;



using rng_type  = mt19937_64;
using rng_int_t = typename rng_type::result_type;
auto constexpr rng_wsize = rng_type::word_size;




auto file_exist( string const& file_name )
    -> bool
{
    auto ifile = ifstream( file_name );
    return ifile.good();
}




// incremental statistic manager.
class statistics
{
private:

    /**
    *******************************************************
    ** member variables.
    *******************************************************
    **/

    double m_mean;
    double m_var;
    double m_max;
    double m_min;
//    double m_median;
    size_t m_count;

public:

    /**
    *******************************************************
    ** constructor.
    *******************************************************
    **/
    statistics() :
        m_mean(0), m_var(0), m_max(0), m_min(0), m_count(0)
    {
        /* no further construction here. */
    }


    /**
    *******************************************************
    ** member functions.
    *******************************************************
    **/
    template < typename T >
    auto update( T const& value )
        -> void
    {
        ++m_count;

        auto prev_mean = m_mean;
        m_mean += (value - prev_mean)/m_count;
        m_var  += (value - prev_mean)*(value - m_mean);

        if( value < m_min ) m_min = value;
        if( value > m_max ) m_max = value;
    }


    template < typename T >
    auto update( vector<T> const& values )
        -> void
    {
        for( auto x : values ) update(x);
    }


    /**
    *******************************************************
    ** friend functions.
    *******************************************************
    **/
    friend
    auto operator<< ( ostream& out, statistics const& stat )
        -> ostream&
    {
        auto sep = '\t';
        out << stat.m_mean << sep
            << stat.m_var/stat.m_count << sep
            << stat.m_max << sep
            << stat.m_min << sep
            << stat.m_max - stat.m_min;

        return out;
    }
};



class cognate
{
public :

    /**
    *******************************************************
    ** typedefs.
    *******************************************************
    **/


private :

    /**
    *******************************************************
    ** member variables.
    *******************************************************
    **/
    /// @dev : map of {word; frequency}.
    unordered_map< size_t, double > m_weight;

public :

    /**
    *******************************************************
    ** constructors.
    *******************************************************
    **/
    cognate() :
        m_weight{ {0,1} } // the common word is '0'.
    {
        /* no further construction here. */
    }

    /**
    *******************************************************
    ** member functions.
    *******************************************************
    **/
    auto sample( rng_type& rng ) const
        -> size_t
    {
        static auto u = uniform_real_distribution< double >( 0, 1 );

        auto uword = u( rng );
        for( auto [word, weight] : m_weight )
            if( uword < weight )
                return word;
            else
                uword -= weight;

        throw runtime_error( "sampling error." );
    }


    auto update_from( cognate const& other,
                      size_t repr,
                      size_t& next_word,
                      size_t max_version,
                      double alpha,
                      double pmut,
                      rng_type& rng )
        -> void
    {
        auto increment = alpha / (1-alpha) / repr;

        // drawing mutation booleans ----------------------
        auto is_mut = vector< int >( other.m_weight.size() );
        auto at_least_one_mutation = false;
        auto bmut = bernoulli_distribution{ pmut };
        for( auto& x : is_mut )
        {
            x = bmut( rng );
            if( x )
                at_least_one_mutation = true;
        }

        if( !at_least_one_mutation )
        {
            for( auto i = size_t{0}; i < repr; ++i )
            {
                auto word = other.sample( rng );
                m_weight[ word ] += increment;
            }
        }
        else
        {
            auto word_correspondance = unordered_map< size_t, size_t >{};
            auto it_mut = is_mut.begin();
            for( auto& x : other.m_weight )
            {
                if( *it_mut )
                    word_correspondance.insert( {x.first, next_word++} );
                else
                    word_correspondance.insert( {x.first, x.first} );

                ++it_mut;
            }

            for( auto i = size_t{0}; i < repr; ++i )
            {
                auto word = word_correspondance[ other.sample( rng ) ];
                m_weight[ word ] += increment;
            }
        }

        // removing the least weights ---------------------
        auto comp = []( auto x, auto y ) { return x.second < y.second; };
        while( m_weight.size() > max_version )
        {
            auto it = min_element( m_weight.begin(), m_weight.end(), comp );
            m_weight.erase( it );
        }

        // weightening the frequencies.
        auto sum = double{0};
        for( auto x : m_weight ) sum += x.second;
        for( auto& x : m_weight ) x.second /= sum;
    }

};


class population
{
public :

    /**
    *******************************************************
    ** typedefs.
    *******************************************************
    **/
    using chr_type = vector< rng_int_t >;

private :

    /**
    *******************************************************
    ** member variables.
    *******************************************************
    **/
    // genetics.
    double const m_gpmut;
    size_t const m_nsnp;
    size_t m_nmale;
    vector< chr_type > m_genetics;

    // linguistics.
    size_t const m_nrl;
    size_t const m_nversion_per_cognate;
    double m_alpha;
    vector< shared_ptr< size_t > > m_next_word;
    vector< size_t > m_repr;
    vector< double > m_lpmut;
    vector< vector < cognate > > m_linguistics;


public :

    /**
    *******************************************************
    ** constructors.
    *******************************************************
    **/
    population(
        size_t n_individual,

        // genetics.
        size_t n_snp,
        double gpmut,

        // linguistics.
        size_t n_cognate,
        size_t n_rl,
        size_t n_version_per_cognate,
        double alpha,
        vector< size_t > repr,
        vector< double > lpmut,

        // rng.
        rng_type& rng
        ) :

        // genetics.
        m_gpmut( gpmut ),
        m_nsnp( n_snp ),
        m_nmale( n_individual / 2 ),
        m_genetics( 2 * n_individual, chr_type( n_snp / rng_wsize + 1 ) ),

        // linguistics.
        m_nrl( n_rl ),
        m_nversion_per_cognate( n_version_per_cognate ),
        m_alpha( alpha ),
        m_next_word( n_cognate ),
        m_repr( move(repr) ),
        m_lpmut( move(lpmut) ),
        m_linguistics( n_individual, vector< cognate >( n_cognate ) )
    {
        // genetics ---------------------------------------
        /// @dev : beta distribution of the frequencies.

        // constructing the truncated distribution probabilities.
        auto a = 4 * n_individual * gpmut;
        auto threshold = double{0.05};
        auto lower     = size_t( n_individual * threshold/2     + 1 );
        auto upper     = size_t( n_individual * (1-threshold/2) );
        auto p = vector< double >( n_individual + 1 - lower - (n_individual-upper) );

        auto cdf = [&] ( size_t x )
            { return boost::math::ibeta( a, a, double(x)/(n_individual+1) ); };

        for( auto i = size_t{0}, end = p.size(); i < end; ++i )
            p[i] = cdf( i + 1 + lower ) - cdf( i + lower );

        auto dist_n1 = discrete_distribution< size_t >( p.begin(), p.end() );

        // assigning.
        auto n1 = vector< size_t >( rng_wsize );
        for( auto ibucket = size_t{0}, end = _n_bucket(); ibucket < end; ++ibucket )
        {
            generate( n1.begin(),
                      n1.end(),
                      [&](){ return 2 * (dist_n1(rng) + lower); } );

            for( auto i = size_t{0}, iend = m_genetics.size(); i < iend; ++i )
            {
                auto& bucket = m_genetics[i][ibucket];

                for( auto bit = size_t{0}; bit < rng_wsize; ++bit )
                    bucket += rng_int_t( i < n1[bit] ) << bit;
            }
        }

        // shuffling.
        /// @dev : the first half individuals are males.
        shuffle( m_genetics.begin(), m_genetics.end(), rng );

        // linguistics ------------------------------------
        for( auto& x : m_next_word ) x = make_shared< size_t >( 1 );
    }


    /**
    *******************************************************
    ** member functions.
    *******************************************************
    **/
    auto event_moran( rng_type& rng )
        -> void
    {
        auto n_individual = _n_individual();

        // drawing indices.
        /// @dev : individual wise index.
        using udist = uniform_int_distribution< size_t >;
        auto i_father = udist{ 0, m_nmale - 1 }( rng );
        auto i_mother = udist{ m_nmale, n_individual - 1 }( rng );
        auto i_offspr = udist{ 0, n_individual - 1 }( rng );

        // genetic inheritance.
        /// @critical : index management : chromosome wise index.
        _genetic_transmission( 2 * i_offspr    , i_father, rng );
        _genetic_transmission( 2 * i_offspr + 1, i_mother, rng );

        // linguistic inheritance.
        _linguistic_transmission( i_offspr, i_father, i_mother, rng );
    }

    auto event_com( rng_type& rng )
        -> void
    {
        auto u = uniform_int_distribution< size_t >{ 0, m_linguistics.size() - 1 };

        // sampling.
        auto a = u( rng ), b = u( rng );
        while( a == b ) b = u( rng );

        auto& listener = m_linguistics[a];
        auto& speaker  = m_linguistics[b];

        for( auto i = size_t{0}, end = _n_cognate(); i < end; ++i )
            listener[i].update_from( speaker[i],
                                     m_repr[i],
                                     *m_next_word[i],
                                     m_nversion_per_cognate,
                                     m_alpha,
                                     m_lpmut[i],
                                     rng );
    }

    // evolve the population during a certain period.
    auto evolve( size_t duration, rng_type& rng )
        -> void
    {
        for( auto t = size_t{0}; t < duration; ++t )
            for( auto i = size_t{0}, end = _n_individual(); i < end; ++i )
            {
                event_moran( rng );

                for( auto j = size_t{0}; j < m_nrl; ++j )
                    event_com( rng );
            }
    }

    // linguistic equilibrium.
    auto wait_L_equilibrium( size_t duration, rng_type& rng )
        -> void
    {
        evolve( duration, rng );
    }


    // demographic events ---------------------------------
    // resize the population (keeping the sexe ratio of 50% - 50%).
    auto resize( size_t n_individual, rng_type& rng )
        -> void
    {
        auto prev_nind = _n_individual();

        if( n_individual != prev_nind )
        {
            // allocating new variables -----------------------
            auto new_nmale = n_individual / 2;
            auto males   = vector< size_t >();
            auto females = vector< size_t >();

            // computing the copy indices ---------------------
            // bottleneck.
            if( n_individual < prev_nind )
            {
                // drawing males.
                /// @dev : individual wise index.
                males.resize( m_nmale );
                iota( males.begin(), males.end(), size_t{0} );
                shuffle( males.begin(), males.end(), rng );
                males.resize( new_nmale );

                // drawing females.
                /// @dev : individual wise index.
                females.resize( prev_nind - m_nmale );
                iota( females.begin(), females.end(), m_nmale );
                shuffle( females.begin(), females.end(), rng );
                females.resize( n_individual - new_nmale );
            }
            // clonal growth.
            else
            {
                using udist = uniform_int_distribution< size_t >;
                auto u_male   = udist{ 0, m_nmale - 1 };
                auto u_female = udist{ m_nmale, prev_nind - 1 };

                males.resize( new_nmale );
                iota( males.begin(), males.begin() + m_nmale, size_t{0} );
                generate( males.begin() + m_nmale      ,
                          males.end()                  ,
                          [&]() { return u_male(rng); }
                          );

                females.resize( n_individual - new_nmale );
                iota( females.begin(), females.begin() + prev_nind - m_nmale, m_nmale );
                generate( females.begin() + prev_nind - m_nmale,
                          females.end()                        ,
                          [&]() { return u_female(rng); }
                          );
            }

            // creating the copy --------------------------
            auto new_gpop = vector< chr_type >( 2 * n_individual );
            auto new_lpop = vector< vector< cognate > >( n_individual );

            auto itg = new_gpop.begin();
            auto itl = new_lpop.begin();

            // "copying" males.
            for( auto m : males )
            {
                *itg++ = ( m_genetics[2*m  ] );
                *itg++ = ( m_genetics[2*m+1] );
                *itl++ = ( m_linguistics[m] );
            }

            // "copying" females.
            for( auto f : females )
            {
                *itg++ = ( m_genetics[2*f  ] );
                *itg++ = ( m_genetics[2*f+1] );
                *itl++ = ( m_linguistics[f] );
            }

            // updating data ------------------------------
            m_genetics.swap( new_gpop );
            m_linguistics.swap( new_lpop );
            m_nmale = new_nmale;
        }
    }




    /**
    *******************************************************
    ** friend functions.
    *******************************************************
    **/
    // genetics.
    friend
    auto G_summary( size_t            n_sim   ,
                    ostream&          out     ,
                    population const& pop     ,
                    size_t            n_sample,
                    rng_type&         rng
                    )
        -> ostream&
    {
        // allocating output variables --------------------
        auto stat_he = statistics{}; // heterogeneity.
        auto n_mono_snp = size_t{0}; // number of monomorphic snps.

        // sampling (indices) -----------------------------
        /// @dev : chromosome wise indices.
        auto u = uniform_int_distribution< size_t >{ 0, pop._n_individual() - 1 };
        auto sample = vector< size_t >( 2 * n_sample );
        for( auto i = size_t{0}; i < n_sample; ++i )
        {
            auto ind = u( rng );
            sample[2*i  ] = 2*ind;
            sample[2*i+1] = 2*ind + 1;
        }


        // computing statistics ---------------------------
        for( auto ibucket = size_t{0}, end = pop._n_bucket(); ibucket < end; ++ibucket )
        {
            // count variable.
            auto count1 = vector< size_t >( rng_wsize );

            // counting 1s.
            for( auto s : sample )
            {
                auto bucket = pop.m_genetics[s][ibucket];
                for( auto bit = size_t{0}; bit < rng_wsize; ++bit )
                    count1[bit] += bucket >> bit & 1;
            }

            // updating statistics.
            /// @dev : taking into account the last bucket.
            auto nlocus = ibucket != end - 1 ? rng_wsize : pop.m_nsnp % rng_wsize;

            for( auto i = size_t{0}; i < nlocus; ++i )
            {
                auto c = count1[i];
                auto f = double(c)/2/n_sample;

                stat_he.update( 1. - f*f - (1.-f)*(1.-f) );

                if( c == 0 || c == 2 * n_sample ) ++n_mono_snp;
            }
        }

        // outputting -------------------------------------
        return out << n_sim << '\t' << n_mono_snp << '\t' << stat_he << '\n';
    }


    friend
    auto G_summary_interpop( size_t            n_sim     ,
                             ostream&          out       ,
                             population const& pop1      ,
                             population const& pop2      ,
                             size_t            n_sample_1,
                             size_t            n_sample_2,
                             rng_type&         rng
                             )
        -> ostream&
    {
        // allocating output variables --------------------
        auto pi_within1 = double{0};
        auto pi_within2 = double{0};
        auto pi_between = double{0};

        // sampling (indices) -----------------------------
        /// @dev : individual wise indices.
        auto u1 = uniform_int_distribution< size_t >{ 0, pop1._n_individual() - 1 };
        auto sample1 = vector< size_t >( n_sample_1 );
        generate_n( sample1.begin(), n_sample_1, [&](){ return u1( rng ); } );

        auto u2 = uniform_int_distribution< size_t >{ 0, pop2._n_individual() - 1 };
        auto sample2 = vector< size_t >( n_sample_2 );
        generate_n( sample2.begin(), n_sample_2, [&](){ return u2( rng ); } );

        // computing statistics ---------------------------
        // first population.
        for( auto ibucket = size_t{0}, end = pop1._n_bucket(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop1.m_nsnp % rng_wsize);

            for( auto i1 = size_t{0}, end1 = sample1.size()-1; i1 < end1; ++i1 )
            {
                auto i1_father_chr = pop1.m_genetics[2*sample1[i1]  ][ibucket];
                auto i1_mother_chr = pop1.m_genetics[2*sample1[i1]+1][ibucket];

                for( auto i2 = i1 + 1, end2 = sample1.size(); i2 < end2; ++i2 )
                {
                    auto i2_father_chr = pop1.m_genetics[2*sample1[i2]  ][ibucket];
                    auto i2_mother_chr = pop1.m_genetics[2*sample1[i2]+1][ibucket];

                    auto pi = __builtin_popcount( (i1_father_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_father_chr ^ i2_mother_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_mother_chr) & mask );

                    pi_within1 += pi;
                    pi_between += pi;
                }
            }
        }

        // second population.
        for( auto ibucket = size_t{0}, end = pop2._n_bucket(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop2.m_nsnp % rng_wsize);

            for( auto i1 = size_t{0}, end1 = sample2.size()-1; i1 < end1; ++i1 )
            {
                auto i1_father_chr = pop1.m_genetics[2*sample2[i1]  ][ibucket];
                auto i1_mother_chr = pop1.m_genetics[2*sample2[i1]+1][ibucket];

                for( auto i2 = i1 + 1, end2 = sample2.size(); i2 < end2; ++i2 )
                {
                    auto i2_father_chr = pop1.m_genetics[2*sample2[i2]  ][ibucket];
                    auto i2_mother_chr = pop1.m_genetics[2*sample2[i2]+1][ibucket];

                    auto pi = __builtin_popcount( (i1_father_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_father_chr ^ i2_mother_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_mother_chr) & mask );

                    pi_within2 += pi;
                    pi_between += pi;
                }
            }
        }

        // only in between.
        //if( pop1._n_individual() != pop2._n_individual() )
          //  throw runtime_error( "non corresponding number of snps." );

        for( auto ibucket = size_t{0}, end = pop1._n_bucket(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop1.m_nsnp % rng_wsize);

            for( auto i1 : sample1 )
            {
                auto i1_father_chr = pop1.m_genetics[2*i1  ][ibucket];
                auto i1_mother_chr = pop1.m_genetics[2*i1+1][ibucket];

                for( auto i2 : sample2 )
                {
                    auto i2_father_chr = pop1.m_genetics[2*i2  ][ibucket];
                    auto i2_mother_chr = pop1.m_genetics[2*i2+1][ibucket];

                    auto pi = __builtin_popcount( (i1_father_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_father_chr) & mask ) +
                              __builtin_popcount( (i1_father_chr ^ i2_mother_chr) & mask ) +
                              __builtin_popcount( (i1_mother_chr ^ i2_mother_chr) & mask );

                    pi_between += pi;
                }
            }
        }

        auto fst = 1. -
            0.5 * ( pi_within1 / n_sample_1 / n_sample_1 + pi_within2 / n_sample_2 / n_sample_2 ) /
            ( pi_between / (n_sample_1 + n_sample_2) / (n_sample_1 + n_sample_2) );

        // outputting -------------------------------------
        return out << n_sim << '\t' << fst << '\n';
    }


    // linguistics.
    friend
    auto L_summary( size_t            n_sim   ,
                    ostream&          out     ,
                    population const& pop     ,
                    size_t            n_sample,
                    rng_type&         rng
                    )
        -> ostream&
    {
        // allocating output variables --------------------
        auto stat_he   = statistics{};
        auto stat_ncog = statistics{}; // number of different cognates.
        auto n_mono_cog = size_t{0};   // number of monomorphic cognates.
//        auto n_uniq_ind = size_t{0};   // number of unique individuals.

        // sampling (indices) -----------------------------
        /// @dev : individual wise indices.
        auto u = uniform_int_distribution< size_t >{ 0, pop._n_individual() - 1 };
        auto sample = vector< size_t >( n_sample );
        generate_n( sample.begin(), n_sample, [&](){ return u( rng ); } );

        // computing statistics ---------------------------
        for( auto icog = size_t{0}, end = pop._n_cognate(); icog < end; ++icog )
        {
            auto count = unordered_map< size_t, size_t >{};

            for( auto s : sample )
            {
                auto c = pop.m_linguistics[s][icog].sample( rng );
                ++count[c];
            }

            // updating.
            stat_ncog.update( count.size() );
            if( count.size() == 1 ) ++n_mono_cog;

            auto ho = double{0};
            for( auto x : count )
            {
                auto f = double(x.second)/n_sample;
                ho += f*f;
            }
            stat_he.update( 1. - ho );
        }

        // outputting -------------------------------------
        return out << n_sim   << '\t' << n_mono_cog << '\t'
                   << stat_he << '\t' << stat_ncog  << '\n';
    }



    friend
    auto L_summary_interpop( size_t            n_sim     ,
                             ostream&          out       ,
                             population const& pop1      ,
                             population const& pop2      ,
                             size_t            n_sample_1,
                             size_t            n_sample_2,
                             rng_type&         rng
                             )
        -> ostream&
    {
        // allocating output variables --------------------
        auto pi_within1 = double{0};
        auto pi_within2 = double{0};
        auto pi_between = double{0};

        // sampling (indices) -----------------------------
        /// @dev : individual wise indices.
        auto u1 = uniform_int_distribution< size_t >{ 0, pop1._n_individual() - 1 };
        auto sample1 = vector< size_t >( n_sample_1 );
        generate_n( sample1.begin(), n_sample_1, [&](){ return u1( rng ); } );

        auto u2 = uniform_int_distribution< size_t >{ 0, pop2._n_individual() - 1 };
        auto sample2 = vector< size_t >( n_sample_2 );
        generate_n( sample2.begin(), n_sample_2, [&](){ return u2( rng ); } );

        // computing statistics ---------------------------
        // first population.
        for( auto icog = size_t{0}, end = pop1._n_cognate(); icog < end; ++icog )
        {
            // first population.
            for( auto i1 = size_t{0}, end1 = sample1.size()-1; i1 < end1; ++i1 )
            {
                for( auto i2 = i1 + 1, end2 = sample1.size(); i2 < end2; ++i2 )
                {
                    auto pi = pop1.m_linguistics[sample1[i1]][icog].sample( rng ) !=
                        pop1.m_linguistics[sample1[i2]][icog].sample( rng );

                    pi_within1 += pi;
                    pi_between += pi;
                }
            }

            // second population.
            for( auto i1 = size_t{0}, end1 = sample2.size()-1; i1 < end1; ++i1 )
            {
                for( auto i2 = i1 + 1, end2 = sample2.size(); i2 < end2; ++i2 )
                {
                    auto pi = pop2.m_linguistics[sample2[i1]][icog].sample( rng ) !=
                        pop2.m_linguistics[sample2[i2]][icog].sample( rng );

                    pi_within2 += pi;
                    pi_between += pi;
                }
            }

            // between populations.
            for( auto i1 : sample1 )
                for( auto i2 : sample2 )
                {
                    auto pi = pop1.m_linguistics[i1][icog].sample( rng ) !=
                        pop2.m_linguistics[i2][icog].sample( rng );

                    pi_between += pi;
                }
        }

        auto fst = 1. -
            0.5 * ( pi_within1 / n_sample_1 / n_sample_1 + pi_within2 / n_sample_2 / n_sample_2 ) /
            ( pi_between / (n_sample_1 + n_sample_2) / (n_sample_1 + n_sample_2) );


        // outputting.
        return out << n_sim << '\t' << fst << '\n';
    }




private :

    /**
    *******************************************************
    ** private functions.
    *******************************************************
    **/
    inline
    auto _n_bucket() const
        -> size_t
    {
        return m_genetics[0].size();
    }

    inline
    auto _n_cognate() const
        -> size_t
    {
        return m_linguistics[0].size();
    }

    inline
    auto _n_individual() const
        -> size_t
    {
        return m_genetics.size()/2;
    }

    auto _genetic_transmission( size_t i_offchr, size_t i_parent, rng_type& rng )
        -> void
    {
        auto& off = m_genetics[i_offchr];
        auto const& father_chr = m_genetics[2*i_parent  ];
        auto const& mother_chr = m_genetics[2*i_parent+1];

        // transmitting.
        /// @dev : loci are independent.
        for( auto i = size_t{0}, end = _n_bucket(); i < end; ++i )
        {
            auto transmit_father_chr = rng();
            off[i] = ( father_chr[i] & transmit_father_chr ) | ( mother_chr[i] & ~transmit_father_chr );
        }

        // muting.
        auto n_mutation = binomial_distribution< size_t >{ m_nsnp, m_gpmut }( rng );
        auto u = uniform_int_distribution< size_t >{ 0, m_nsnp - 1 };

        for( auto i = size_t{0}; i < n_mutation; ++i )
        {
            // drawing the mutation position.
            auto mut_pos = u(rng);

            // getting the bucket and internal position.
            auto bucket = mut_pos / rng_wsize;
            auto ipos = mut_pos % rng_wsize;

            // changing the adequate bit.
            off[bucket] ^= (rng_int_t{1} << ipos);
        }
    }


    // linguistic inheritance.
    auto _linguistic_transmission( size_t i_offspr, size_t i_father, size_t i_mother, rng_type& rng )
        -> void
    {
        // one of the parent is the teacher.
        auto i_teacher = rng() & 1 ? i_father : i_mother;

        m_linguistics[i_offspr] = m_linguistics[i_teacher];
    }
};





/*
 * A class to read data from a csv file.
 */
class CSVReader
{
private :

    string m_filename;
    string m_delimiter;

public:

    CSVReader( string filename, string delm = "," ) :
        m_filename( filename ), m_delimiter( delm )
    {
        /* no further construction here. */
    }

    // function to fetch data from a CSV File
    auto getData() const
        -> vector< vector< string > >
    {
        auto dataList = vector< vector< string > >{};

        if( auto file = ifstream( m_filename ); file.fail() )
            runtime_error( "input error." );
        else
        {
            auto line = string{};
            while( getline(file, line) )
            {
                auto splitted_string = vector< string >{};
                boost::algorithm::split( splitted_string, line, boost::is_any_of( m_delimiter ) );
                dataList.push_back( splitted_string );
            }
        }

        return dataList;
    }

};




/**********************************************************
** main
**********************************************************/
int main( int argc, char** argv )
{
    try
    {
        /**************************************************
        ** options ****************************************
        **************************************************/
        // random number generator ------------------------
        auto rng = rng_type{ random_device{}() };
//        auto rng = rng_type{ 0 };

        // linguistics ------------------------------------
        auto n_cognate = size_t{34};
        auto n_version_per_cognate = size_t{4};

        // genetics ---------------------------------------
        auto n_snp    = size_t{85425};

        // simulation parameters --------------------------
        auto n_sample_0 = size_t{23};
        auto n_sample_1 = size_t{10};
        auto g_time     = size_t{25} ;

        auto n_sim      = size_t{0} ;

        auto n_Rg     = size_t{1} ;
        auto n_Rl     = size_t{10} ;
        auto N0       = size_t{100} ;
        auto N1_ST    = size_t{50} ;
        auto N2_ST    = size_t{100} ;
        auto N1_MA    = size_t{100} ;
        auto N2_MA    = size_t{100} ;
        auto N1_ST_MA = size_t{200} ;
        auto N2_ST_MA = size_t{200} ;
        auto Nb_MA    = size_t{100} ;
        auto t0       = size_t{300} ;
        auto t_MA_b   = size_t{200} ;
        auto t_MA_e   = size_t{100} ;
        auto t_Nexp   = size_t{50} ;
        auto L_repr_min = size_t{1} ;
        auto L_repr_max = size_t{1} ;

        auto L_alpha  = double{0.5} ;
        auto L_mut    = double{0.001} ;
        auto G_mut    = double{0.00000001} ;

        // parsing.
        if( argc > 1 )
        {
            n_sim = stoul( argv[1] );

            // string representation of variables.
            auto string_repr_ul = map< string, size_t* >{
                { "n_Rg"      , &n_Rg       },
                { "n_Rl"      , &n_Rl       },
                { "N0"        , &N0         },
                { "N1_ST"     , &N1_ST      },
                { "N2_ST"     , &N2_ST      },
                { "N1_MA"     , &N1_MA      },
                { "N2_MA"     , &N2_MA      },
                { "N1_ST_MA"  , &N1_ST_MA   },
                { "N2_ST_MA"  , &N2_ST_MA   },
                { "Nb_MA"     , &Nb_MA      },
                { "t0"        , &t0         },
                { "t_MA_b"    , &t_MA_b     },
                { "t_MA_e"    , &t_MA_e     },
                { "t_Nexp"    , &t_Nexp     },
                { "L_repr_min", &L_repr_min },
                { "L_repr_max", &L_repr_max },
            };

            auto string_repr_double = map< string, double* >{
                { "L_alpha", &L_alpha },
                { "L_mut"  , &L_mut   },
                { "G_mut"  , &G_mut   }
            };


            // Creating an object of CSVWriter
            auto reader = CSVReader( "param.txt", "\t" );

            // Get the data from CSV File
            auto dataList = reader.getData();

            auto param_label = dataList[0];
            auto param_value = dataList[ stoul( argv[1] ) ];

            for( auto i = size_t{0}, end = param_label.size(); i < end ; ++i )
                if( string_repr_ul.find( param_label[i]) != string_repr_ul.end() )
                    *string_repr_ul[ param_label[i] ] = stoul( param_value[i] );
                else if( string_repr_double.find( param_label[i]) != string_repr_double.end() )
                    *string_repr_double[ param_label[i] ] = stod( param_value[i] );
                else
                    throw runtime_error( param_label[i] + " not provided in parameter file." );
        }


        /**************************************************
        ** priors *****************************************
        **************************************************/
        // linguistics ------------------------------------
        // representability
        /// @info : importance of the word in a speech.
        auto u_repr = uniform_int_distribution< size_t >{ L_repr_min, L_repr_max };
        auto L_repr = vector< size_t >( n_cognate );
        for( auto& x : L_repr ) x = u_repr( rng );

        // mutation rates
        auto L_pmut = vector< double >( n_cognate );

        auto alpha = double{2};
        using gdist = gamma_distribution< double >;
        auto gammaAlphaL = gdist{ alpha, 1. };
        auto gammaBetaL  = gdist{ alpha * (1. - L_mut)/L_mut, 1. };

        for( auto& x : L_pmut )
        {
            auto g = gammaAlphaL( rng );
            x = g / ( g + gammaBetaL( rng ) );
        }

        // genetics ---------------------------------------
        auto G_pmut = G_mut;

        // Checking parameters.
        assert( L_repr.size() == n_cognate );
        assert( L_pmut.size() == n_cognate );
        assert( n_sample_0 < N2_ST_MA );

        /**************************************************
        ** simulation *************************************
        **************************************************/
        // initialising a population ----------------------
        auto pop = population( N0,
                               // genetics.
                               n_snp,
                               G_pmut,

                               // linguistics.
                               n_cognate,
                               n_Rl,
                               n_version_per_cognate,
                               L_alpha,
                               L_repr,
                               L_pmut,
                               rng
                               );

        // scenario ---------------------------------------
        // waiting for equilibra.
        pop.wait_L_equilibrium( 10 * N0 / n_Rl / L_alpha, rng );

        // expansion 1 ST
        pop.resize( N1_ST, rng );

        // evolution ST
        pop.evolve( (t0 - t_MA_b) / g_time, rng );

        // split ST / MA
        auto pop1 = pop ;

        //////////////////////
        /// Evolution MA /////

        // bottleneck MA
        pop1.resize( Nb_MA, rng );

        // evolution MA
        pop1.evolve( (t_MA_b - t_MA_e) / g_time, rng );

        // end bottleneck MA
        pop1.resize( N1_MA, rng );

        // evolution MA
        pop1.evolve( (t_MA_e - t_Nexp) / g_time, rng );

        // expansion MA
        pop1.resize( N2_MA, rng );

        // evolution MA
        pop1.evolve( t_Nexp / g_time, rng );


        ////////////////////////
        ///// Evolution ST /////

        // evolution ST
        pop.evolve( (t_MA_b - t_Nexp) / g_time, rng );

        // expansion 2 ST
        pop.resize( N2_ST, rng );

        // evolution ST
        pop.evolve( t_Nexp / g_time, rng );



        /**************************************************
        ** outputting *************************************
        **************************************************/
        auto header_G = "n_sim\tS\tmeanD\tvarD\tmaxD\tminD\trangeD\n" ;
        auto header_L = "n_sim\tS\tmeanD\tvarD\tmaxD\tminD\trangeD\t"
                         "meanR\tvarR\tmaxR\tminR\trangeR\n" ;

        auto G_filename_0 = "stat_G_out_0" ;
        auto G_file_existed_0 = file_exist(G_filename_0);

        if( auto G_out = ofstream( G_filename_0, ios::app ); G_out.fail() )
            throw runtime_error( "cannot open summary file." );
        else
        {
            if( !G_file_existed_0 )
                G_out << header_G;

            G_summary( n_sim, G_out, pop, n_sample_0, rng );
        }

        auto L_filename_0 = "stat_L_out_0";
        auto L_file_existed_0 = file_exist( L_filename_0 );
        if( auto L_out = ofstream( L_filename_0, ios::app ); L_out.fail() )
            throw runtime_error( "cannot open summary file." );
        else
        {
            if( !L_file_existed_0 )
                L_out << header_L;

            L_summary( n_sim, L_out, pop, n_sample_0, rng );
        }


        auto G_filename_1 = "stat_G_out_1" ;
        auto G_file_existed_1 = file_exist( G_filename_1 );
        if( auto G_out = ofstream( G_filename_1, ios::app ); G_out.fail() )
            throw runtime_error( "cannot open summary file." );
        else
        {
            if( !G_file_existed_1 )
                G_out << header_G;

            G_summary( n_sim, G_out, pop1, n_sample_1, rng );
        }

        auto L_filename_1 = "stat_L_out_1";
        auto L_file_existed_1 = file_exist( L_filename_1 );
        if( auto L_out = ofstream( L_filename_1, ios::app ); L_out.fail() )
            throw runtime_error( "cannot open summary file." );
        else
        {
            if( !L_file_existed_1 )
                L_out << header_L;

            L_summary( n_sim, L_out, pop1, n_sample_1, rng );
        }

        // writing interpopulation summaries. /////////////
        auto header_pair = "n_sim\tFst\n";

        auto G_filename_0_1 = "stat_G_out_0_1";
            auto G_file_existed_0_1 = file_exist( G_filename_0_1 );
            if( auto G_out = ofstream( G_filename_0_1, ios::app ); G_out.fail() )
                throw runtime_error( "cannot open summary file." );
            else
            {
                if( !G_file_existed_0_1 )
                    G_out << header_pair;

                G_summary_interpop( n_sim,
                                    G_out,
                                    pop,
                                    pop1,
                                    n_sample_0,
                                    n_sample_1,
                                    rng );
            }

            auto L_filename_0_1 = "stat_L_out_0_1";
            auto L_file_existed_0_1 = file_exist( L_filename_0_1 );
            if( auto L_out = ofstream( L_filename_0_1, ios::app ); L_out.fail() )
                throw runtime_error( "cannot open summary file." );
            else
            {
                if( !L_file_existed_0_1 )
                    L_out << header_pair;

                L_summary_interpop( n_sim,
                                    L_out,
                                    pop,
                                    pop1,
                                    n_sample_0,
                                    n_sample_1,
                                    rng );
            }

    }
    /******************************************************
    ** catching exceptions ********************************
    ******************************************************/
    catch( std::exception& e )
    {
        cerr << "plgs: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
