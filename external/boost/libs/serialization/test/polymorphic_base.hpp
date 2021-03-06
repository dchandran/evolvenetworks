#ifndef POLYMORPHIC_BASE_HPP
#define POLYMORPHIC_BASE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// polymorphic_base.hpp    simple class test

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/serialization/access.hpp>
#include <boost/serialization/type_info_implementation.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>

class polymorphic_base
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(
        Archive & /* ar */, 
        const unsigned int /* file_version */
    ){}
public:
    // note that since this class uses the "no_rtti"
    // extended_type_info implementation, it MUST
    // implement this function
    virtual const char * get_key() const;
    virtual ~polymorphic_base(){};
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(polymorphic_base)

BOOST_CLASS_TYPE_INFO(
    polymorphic_base,
    boost::serialization::extended_type_info_no_rtti<polymorphic_base>
)

#endif // POLYMORPHIC_BASE_HPP
