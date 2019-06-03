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
    unordered_map< size_t, double > m_freq;


public :

    /**
    *******************************************************
    ** constructors.
    *******************************************************
    **/
    cognate() :
        m_freq{ {0,1} } // the common word is '0'.
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
        for( auto& x : m_freq )
            if( uword < x.second )
                return x.first;
            else
                uword -= x.second;

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
        // weighting the listener.
        for( auto& x : m_freq ) x.second *= 1. - alpha;

        auto mut = vector< int >( other.m_freq.size() );
        auto is_mut = false;
        auto bmut = bernoulli_distribution( pmut );
        for( auto& x : mut )
        {
            x = bmut( rng );
            if( x == 1 )
                is_mut = true;
        }

        if( is_mut )
        {
            auto word_correspondance = unordered_map< size_t, size_t >{};
            auto it_mut = mut.begin();
            for( auto& x : other.m_freq )
            {
                if( *it_mut )
                    word_correspondance.insert( {x.first, next_word++} );
                else
                    word_correspondance.insert( {x.first, x.first});

                ++it_mut;
            }

            for( auto i = size_t{0}; i < repr; ++i )
            {
                auto word = word_correspondance[ other.sample( rng ) ];

                if( m_freq.find( word ) != m_freq.end() )
                    m_freq[ word ] += alpha / repr;
                else
                    m_freq.insert( {word, alpha / repr} );
            }
        }
        else
        {
            for( auto i = size_t{0}; i < repr; ++i )
            {
                auto word = other.sample( rng );

                if( m_freq.find( word ) != m_freq.end() )
                    m_freq[ word ] += alpha / repr;
                else
                    m_freq.insert( {word, alpha / repr} );
            }
        }

        // removing the least weights ---------------------
        auto comp = []( auto x, auto y ) -> bool { return x.second < y.second; };
        while( m_freq.size() > max_version )
        {
            auto it = min_element( m_freq.begin(), m_freq.end(), comp );
            m_freq.erase( it );
        }

        // weightening the frequencies.
        auto sum = double{0};
        for( auto x : m_freq ) sum += x.second;
        for( auto& x : m_freq ) x.second /= sum;
    }


private :

    /**
    *******************************************************
    ** private functions.
    *******************************************************
    **/

};




class population
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

    size_t m_nindividual;

    // genetics.
    /// @dev : index management is individual wise at
    ///        first, then is transformed appropriately to
    ///        take account of the diploidy.
    size_t m_nmale; // excluding end index of male.
    size_t m_nsnp;
    double m_Gpmut;

    /// @opt : the internal vector is a pointer, thus the
    ///        size is not managed automatically anymore.
    ///        A good index management is hence critical.
    /// @dev : diploids.
    /// @dev : first vector is for loci,; the second for
    ///        individuals.
    vector< unique_ptr< rng_int_t[] > > m_genetics;

    // linguistics.
    size_t m_nrl;
    size_t m_ncognate;
    size_t m_nversion_per_cognate;
    vector< shared_ptr< size_t > > m_next_word;
    vector< size_t > m_repr;
    double m_alpha;
    vector< double > m_Lpmut;

    /// @dev : first vector is for cognates; the second for
    ///        individuals.
    vector< unique_ptr< cognate[] > > m_linguistics;


public :

    /**
    *******************************************************
    ** constructors.
    *******************************************************
    **/
    // main constructor.
    population( size_t n_individual, size_t n_snp, double G_pmut,
                size_t n_rl, size_t n_cognate, size_t n_vpercognate,
                vector< size_t > const& repr, double alpha, vector< double > const& L_pmut ) :
        m_nindividual( n_individual ),

        // genetics.
        m_nmale( n_individual/2 ),
        m_nsnp( n_snp ),
        m_Gpmut( G_pmut ),
        m_genetics( n_snp / rng_wsize + 1 ),

        // linguistics.
        m_nrl( n_rl ),
        m_ncognate( n_cognate ),
        m_nversion_per_cognate( n_vpercognate ),
        m_next_word( n_cognate ),
        m_repr( repr ),
        m_alpha( alpha ),
        m_Lpmut( L_pmut ),
        m_linguistics( n_cognate )
    {
        // genetics.
        for( auto& x : m_genetics )
            x = make_unique< rng_int_t[] >( 2 * n_individual ); // diploids.

        // linguistics.
        for( auto& x : m_next_word )
            x = make_shared< size_t >( 2 );

        for( auto& x : m_linguistics )
            x = make_unique< cognate[] >( n_individual );
    }

    // copy constructor.
    population( population const& other ) :
        m_nindividual( other.m_nindividual ),

        // genetics.
        m_nmale( other.m_nmale ),
        m_nsnp( other.m_nsnp ),
        m_Gpmut( other.m_Gpmut ),
        m_genetics( other.m_nsnp / rng_wsize + 1 ),

        // linguistics.
        m_nrl( other.m_nrl ),
        m_ncognate( other.m_ncognate ),
        m_nversion_per_cognate( other.m_nversion_per_cognate ),
        m_next_word( other.m_next_word ),
        m_repr( other.m_repr ),
        m_alpha( other.m_alpha ),
        m_Lpmut( other.m_Lpmut ),
        m_linguistics( other.m_ncognate )
    {
        // genetics.
        for( auto i = size_t{0}, end = m_genetics.size(); i < end; ++i )
        {
            m_genetics[i] = make_unique< rng_int_t[] >( 2 * m_nindividual );
            for( auto j = size_t{0}; j < 2 * m_nindividual; ++j )
                m_genetics[i][j] = other.m_genetics[i][j];
        }

        // linguistics.
        for( auto i = size_t{0}, end = m_linguistics.size(); i < end; ++i )
        {
            m_linguistics[i] = make_unique< cognate[] >( m_nindividual );
            for( auto j = size_t{0}; j < m_nindividual; ++j )
                m_linguistics[i][j] = other.m_linguistics[i][j];
        }
    }


    /**
    *******************************************************
    ** member functions.
    *******************************************************
    **/
    // unit events ----------------------------------------
    auto event_moran( rng_type& rng )
        -> void
    {
        // drawing indices.
        /// @dev : individual wise index.
        using udist = uniform_int_distribution< size_t >;
        auto i_father = udist{ 0, m_nmale - 1 }( rng );
        auto i_mother = udist{ m_nmale, m_nindividual - 1 }( rng );
        auto i_offspr = udist{ 0, m_nindividual - 1 }( rng );

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
        auto u = uniform_int_distribution< size_t >{ 0, m_nindividual - 1 };

        // sampling.
        auto a = u( rng ), b = u( rng );
        while( a == b ) b = u( rng );

        for( auto i = size_t{0}; i < m_ncognate; ++i )
            m_linguistics[i][a].update_from( m_linguistics[i][b],
                                             m_repr[i],
                                             *m_next_word[i],
                                             m_nversion_per_cognate,
                                             m_alpha,
                                             m_Lpmut[i],
                                             rng );
    }


    // demographic events ---------------------------------
    // resize the population (keeping the sexe ratio of 50% - 50%).
    auto resize( size_t n_individual, rng_type& rng )
        -> void
    {
        if( n_individual != m_nindividual )
        {
            // allocating new variables -----------------------
            auto new_nmale = n_individual / 2;
            auto males   = vector< size_t >();
            auto females = vector< size_t >();

            // computing the copy indices ---------------------
            // bottleneck.
            if( n_individual < m_nindividual )
            {
                // drawing males.
                /// @dev : individual wise index.
                males.resize( m_nmale );
                iota( males.begin(), males.end(), size_t{0} );
                shuffle( males.begin(), males.end(), rng );

                // drawing females.
                /// @dev : individual wise index.
                females.resize( m_nindividual - m_nmale );
                iota( females.begin(), females.end(), m_nmale );
                shuffle( females.begin(), females.end(), rng );
            }
            // clonal growth.
            else if( n_individual > m_nindividual )
            {
                using udist = uniform_int_distribution< size_t >;
                auto u_male   = udist{ 0, m_nmale - 1 };
                auto u_female = udist{ m_nmale, m_nindividual - 1 };

                males.resize( new_nmale );
                iota( males.begin(), males.begin() + m_nmale, size_t{0} );
                generate( males.begin() + m_nmale      ,
                          males.end()                  ,
                          [&]() { return u_male(rng); }
                          );

                females.resize( n_individual - new_nmale );
                iota( females.begin(), females.begin() + m_nindividual - m_nmale, m_nmale );
                generate( females.begin() + m_nindividual - m_nmale,
                          females.end()                            ,
                          [&]() { return u_female(rng); }
                          );
            }

            // updating genetic variables ---------------------
            for( auto& bucket : m_genetics )
            {
                // allocating the new bucket.
                auto new_bucket = make_unique< rng_int_t[] >( 2 * n_individual );
                auto it = new_bucket.get();

                // copying males.
                for( auto i = size_t{0}; i < new_nmale; ++i )
                {
                    /// @critical : index management : chromosome wise index.
                    *it++ = bucket[ 2*males[i]     ];
                    *it++ = bucket[ 2*males[i] + 1 ];
                }

                // copying females.
                for( auto i = size_t{0}; i < n_individual - new_nmale; ++i )
                {
                    /// @critical : index management : chromosome wise index.
                    *it++ = bucket[ 2*females[i]     ];
                    *it++ = bucket[ 2*females[i] + 1 ];
                }

                // swapping buckets.
                bucket.swap( new_bucket );
            }

            // updating linguistic variables ------------------
            for( auto& cog : m_linguistics )
            {
                // allocating the new bucket.
                auto new_cog = make_unique< cognate[] >( n_individual );
                auto it = new_cog.get();

                // copying males.
                /// @critical : index management.
                for( auto i = size_t{0}; i < new_nmale; ++i )
                    *it++ = cog[ males[i] ];

                // copying females.
                /// @critical : index management.
                for( auto i = size_t{0}; i < n_individual - new_nmale; ++i )
                    *it++ = cog[ females[i] ];

                // swapping buckets.
                cog.swap( new_cog );
            }

            // updating size related variables ----------------
            m_nindividual = n_individual;
            m_nmale = new_nmale;
        }
    }


    // evolve the population during a certain period.
    auto evolve( size_t duration, rng_type& rng )
        -> void
    {
        for( auto t = size_t{0}; t < duration; ++t )
            for( auto i = size_t{0}; i < m_nindividual; ++i )
            {
                event_moran( rng );

                for( auto j = size_t{0}; j < m_nrl; ++j )
                    event_com( rng );
            }
    }


    // genetic equilibrium
    /// @todo : incorporate the duration?
    auto wait_G_equilibrium( size_t duration, rng_type& rng )
        -> void
    {
        for( auto i = size_t{0}; i < duration; ++i )
            for( auto j = size_t{0}; j < m_nindividual; ++j )
                event_moran( rng );
    }


    // linguistic equilibrium.
    /// @todo : incorporate the duration?
    auto wait_L_equilibrium( size_t duration, rng_type& rng )
        -> void
    {
        evolve( duration, rng );
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
        /// @dev : individual wise indices.
        auto u = uniform_int_distribution< size_t >{ 0, pop.m_nindividual - 1 };
        auto sample = vector< size_t >( n_sample );
        generate_n( sample.begin(), n_sample, [&](){ return u( rng ); } );

        // computing statistics ---------------------------
        for( auto ibucket = size_t{0}, end = pop.m_genetics.size(); ibucket < end; ++ibucket )
        {
            // count variable.
            auto& bucket = pop.m_genetics[ibucket];
            auto count1 = vector< size_t >( rng_wsize );

            // counting 1s.
            for( auto s : sample )
            {
                /// @dev : chomosome wise index.
                auto father_chr = bucket[2*s  ];
                auto mother_chr = bucket[2*s+1];

                for( auto bit = size_t{0}; bit < rng_wsize; ++bit )
                    count1[bit] += ( father_chr >> bit & 1 ) + ( mother_chr >> bit & 1 );
            }


            // updating statistics.
            /// @dev : taking into account the last bucket.
            auto nlocus = ibucket != end - 1? rng_wsize : pop.m_nsnp % rng_wsize;

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
        auto u1 = uniform_int_distribution< size_t >{ 0, pop1.m_nindividual - 1 };
        auto sample1 = vector< size_t >( n_sample_1 );
        generate_n( sample1.begin(), n_sample_1, [&](){ return u1( rng ); } );

        auto u2 = uniform_int_distribution< size_t >{ 0, pop2.m_nindividual - 1 };
        auto sample2 = vector< size_t >( n_sample_2 );
        generate_n( sample2.begin(), n_sample_2, [&](){ return u2( rng ); } );

        // computing statistics ---------------------------
        // first population.
        for( auto ibucket = size_t{0}, end = pop1.m_genetics.size(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop1.m_nsnp % rng_wsize);

            for( auto i1 = size_t{0}, end1 = sample1.size()-1; i1 < end1; ++i1 )
            {
                auto i1_father_chr = pop1.m_genetics[ibucket][2*sample1[i1]  ];
                auto i1_mother_chr = pop1.m_genetics[ibucket][2*sample1[i1]+1];

                for( auto i2 = i1 + 1, end2 = sample1.size(); i2 < end2; ++i2 )
                {
                    auto i2_father_chr = pop1.m_genetics[ibucket][2*sample1[i2]  ];
                    auto i2_mother_chr = pop1.m_genetics[ibucket][2*sample1[i2]+1];

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
        for( auto ibucket = size_t{0}, end = pop2.m_genetics.size(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop2.m_nsnp % rng_wsize);

            for( auto i1 = size_t{0}, end1 = sample2.size()-1; i1 < end1; ++i1 )
            {
                auto i1_father_chr = pop1.m_genetics[ibucket][2*sample2[i1]  ];
                auto i1_mother_chr = pop1.m_genetics[ibucket][2*sample2[i1]+1];

                for( auto i2 = i1 + 1, end2 = sample2.size(); i2 < end2; ++i2 )
                {
                    auto i2_father_chr = pop1.m_genetics[ibucket][2*sample2[i2]  ];
                    auto i2_mother_chr = pop1.m_genetics[ibucket][2*sample2[i2]+1];

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
        if( pop1.m_genetics.size() != pop2.m_genetics.size() )
            throw runtime_error( "non corresponding number of snps." );

        for( auto ibucket = size_t{0}, end = pop1.m_genetics.size(); ibucket < end; ++ibucket )
        {
            auto mask = size_t(-1);
            if( ibucket == end - 1 ) // masking for least bits.
                mask = ~(mask << pop1.m_nsnp % rng_wsize);

            for( auto i1 : sample1 )
            {
                auto i1_father_chr = pop1.m_genetics[ibucket][2*i1  ];
                auto i1_mother_chr = pop1.m_genetics[ibucket][2*i1+1];

                for( auto i2 : sample2 )
                {
                    auto i2_father_chr = pop1.m_genetics[ibucket][2*i2  ];
                    auto i2_mother_chr = pop1.m_genetics[ibucket][2*i2+1];

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
        auto u = uniform_int_distribution< size_t >{ 0, pop.m_nindividual - 1 };
        auto sample = vector< size_t >( n_sample );
        generate_n( sample.begin(), n_sample, [&](){ return u( rng ); } );

        // computing statistics ---------------------------
        for( auto& cog : pop.m_linguistics )
        {
            auto count = unordered_map< size_t, size_t >{};

            for( auto s : sample )
            {
                auto c = cog[s].sample( rng );
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
        auto u1 = uniform_int_distribution< size_t >{ 0, pop1.m_nindividual - 1 };
        auto sample1 = vector< size_t >( n_sample_1 );
        generate_n( sample1.begin(), n_sample_1, [&](){ return u1( rng ); } );

        auto u2 = uniform_int_distribution< size_t >{ 0, pop2.m_nindividual - 1 };
        auto sample2 = vector< size_t >( n_sample_2 );
        generate_n( sample2.begin(), n_sample_2, [&](){ return u2( rng ); } );

        // computing statistics ---------------------------
        // first population.
        for( auto& cog : pop1.m_linguistics )
        {
            for( auto i1 = size_t{0}, end1 = sample1.size()-1; i1 < end1; ++i1 )
            {
                for( auto i2 = i1 + 1, end2 = sample1.size(); i2 < end2; ++i2 )
                {
                    auto pi = cog[sample1[i1]].sample( rng ) != cog[sample1[i2]].sample( rng );

                    pi_within1 += pi;
                    pi_between += pi;
                }
            }
        }

        // second population.
        for( auto& cog : pop2.m_linguistics )
        {
            for( auto i1 = size_t{0}, end1 = sample2.size()-1; i1 < end1; ++i1 )
            {
                for( auto i2 = i1 + 1, end2 = sample2.size(); i2 < end2; ++i2 )
                {
                    auto pi = cog[sample2[i1]].sample( rng ) != cog[sample2[i2]].sample( rng );

                    pi_within2 += pi;
                    pi_between += pi;
                }
            }
        }

        // only in between.
        if( pop1.m_linguistics.size() != pop1.m_linguistics.size() )
            throw runtime_error( "not matching number of cognates." );

        for( auto i = size_t{0}, end = pop1.m_linguistics.size(); i < end; ++i )
        {
            for( auto i1 : sample1 )
                for( auto i2 : sample2 )
                {
                    auto pi = pop1.m_linguistics[i][i1].sample( rng ) !=
                        pop2.m_linguistics[i][i2].sample( rng );

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
    // genetic inheritance.
    /// @dev : %i_offchr is a chromosome wise index.
    /// @dev : %i_parent is a individual wise index.
    auto _genetic_transmission( size_t i_offchr, size_t i_parent, rng_type& rng )
        -> void
    {
        // transmitting.
        /// @dev : loci are independent.
        for( auto ibucket = size_t{0}, end = m_genetics.size(); ibucket < end; ++ibucket )
        {
            auto& bucket = m_genetics[ibucket];

            auto transmit_fat_chr = rng();
            bucket[i_offchr] = ( bucket[2*i_parent  ] & transmit_fat_chr ) |
                               ( bucket[2*i_parent+1] & ~transmit_fat_chr );
        }

        // muting.
        auto n_mutation = binomial_distribution< size_t >{ m_nsnp, m_Gpmut }( rng );
        auto u = uniform_int_distribution< size_t >{ 0, m_nsnp - 1 };

        for( auto i = size_t{0}; i < n_mutation; ++i )
        {
            // drawing the mutation position.
            auto mut_pos = u(rng);

            // getting the bucket and internal position.
            auto bucket = mut_pos / rng_wsize;
            auto ipos = mut_pos % rng_wsize;

            // changing the adequate bit.
            m_genetics[bucket][i_offchr] ^= (rng_int_t{1} << ipos);
        }
    }

    // linguistic inheritance.
    auto _linguistic_transmission( size_t i_offspr, size_t i_father, size_t i_mother, rng_type& rng )
        -> void
    {
        // one of the parent is the teacher.
        auto i_teacher = rng() & 1 ? i_father : i_mother;

        for( auto& cog : m_linguistics )
            cog[i_offspr] = cog[i_teacher];
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

        // linguistics ------------------------------------
        auto n_cognate = size_t{34};
        auto n_version_per_cognate = size_t{4};

        // genetics ---------------------------------------
        auto n_snp    = size_t{100};

        // simulation parameters --------------------------
        auto n_sample_0 = size_t{24};
        auto n_sample_1 = size_t{10};
        
        auto g_time     = size_t{25} ;

        auto n_sim      = size_t{0} ;

        auto n_Rg     = size_t{1} ;
        auto n_Rl     = size_t{10} ;
        auto N0       = size_t{100} ;
        auto N1_ST    = size_t{100} ;
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
        auto G_mut    = double{0.00001} ;

        // parsing.
        if( argc > 1 )
        {
            n_snp = size_t{2};
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
        auto gammaAlphaG = gdist{ alpha, 1. };
        auto gammaBetaG  = gdist{ alpha*(1 - G_mut) / G_mut, 1. };
        auto g = gammaAlphaG( rng );

        auto G_pmut = g / ( g + gammaBetaG( rng ) );

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
                               n_Rl,
                               n_cognate,
                               n_version_per_cognate,
                               L_repr,
                               L_alpha,
                               L_pmut
                               );

        // scenario ---------------------------------------
        // waiting for equilibra.
        pop.wait_G_equilibrium( 12 * N0, rng );
        pop.wait_L_equilibrium( 10 * N0 / n_Rl / L_alpha, rng );


        // expansion 1 ST
        pop.resize( N1_ST_MA, rng );

        // evolution ST
        pop.evolve( (t0 - t_Nexp) / g_time, rng );

        // expansion 2 ST
        pop.resize( N2_ST_MA, rng );

        // evolution ST
        pop.evolve( t_Nexp/g_time, rng );



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

            G_summary( n_sim, G_out, pop, n_sample_1, rng );
        }

        auto L_filename_1 = "stat_L_out_1";
        auto L_file_existed_1 = file_exist( L_filename_1 );
        if( auto L_out = ofstream( L_filename_1, ios::app ); L_out.fail() )
            throw runtime_error( "cannot open summary file." );
        else
        {
            if( !L_file_existed_1 )
                L_out << header_L;

            L_summary( n_sim, L_out, pop, n_sample_1, rng );
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
                                    pop,
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
                                    pop,
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
