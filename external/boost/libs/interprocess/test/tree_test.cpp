//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2004-2007. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/interprocess for documentation.
//
//////////////////////////////////////////////////////////////////////////////
#include <boost/interprocess/detail/config_begin.hpp>
#include <set>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/set.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/indexes/map_index.hpp>
#include <boost/interprocess/indexes/iset_index.hpp>
#include <boost/interprocess/mem_algo/simple_seq_fit.hpp>
#include "print_container.hpp"
#include "movable_int.hpp"
#include "dummy_test_allocator.hpp"
#include "set_test.hpp"
#include "map_test.hpp"
#include "emplace_test.hpp"

///////////////////////////////////////////////////////////////////
//                                                               //
//  This example repeats the same operations with std::set and   //
//  shmem_set using the node allocator                           //
//  and compares the values of both containers                   //
//                                                               //
///////////////////////////////////////////////////////////////////

using namespace boost::interprocess;

//Customize managed_shared_memory class
typedef basic_managed_shared_memory
   <char,
    simple_seq_fit<mutex_family, offset_ptr<void> >,
    map_index
   > my_managed_shared_memory;

//We will work with narrow characters for shared memory objects
//Alias an integer node allocator type
typedef allocator<int, my_managed_shared_memory::segment_manager>   
   shmem_allocator_t;
typedef allocator<std::pair<const int, int>, my_managed_shared_memory::segment_manager>   
   shmem_node_pair_allocator_t;
typedef allocator<test::movable_int, my_managed_shared_memory::segment_manager>   
   shmem_movable_allocator_t;
typedef allocator<std::pair<const test::movable_int, test::movable_int>, my_managed_shared_memory::segment_manager>   
   shmem_movable_node_pair_allocator_t;
typedef allocator<test::movable_and_copyable_int, my_managed_shared_memory::segment_manager>   
   shmem_move_copy_allocator_t;
typedef allocator<std::pair<const test::movable_and_copyable_int, test::movable_and_copyable_int>, my_managed_shared_memory::segment_manager>   
   shmem_move_copy_node_pair_allocator_t;

//Alias standard types
typedef std::set<int>                                          MyStdSet;
typedef std::multiset<int>                                     MyStdMultiSet;
typedef std::map<int, int>                                     MyStdMap;
typedef std::multimap<int, int>                                MyStdMultiMap;

//Alias non-movable types
typedef set<int, std::less<int>, shmem_allocator_t>       MyShmSet;
typedef multiset<int, std::less<int>, shmem_allocator_t>  MyShmMultiSet;
typedef map<int, int, std::less<int>, shmem_node_pair_allocator_t>  MyShmMap;
typedef multimap<int, int, std::less<int>, shmem_node_pair_allocator_t>  MyShmMultiMap;

//Alias movable types
typedef set<test::movable_int, std::less<test::movable_int>
            ,shmem_movable_allocator_t>                   MyMovableShmSet;
typedef multiset<test::movable_int, 
      std::less<test::movable_int>, 
      shmem_movable_allocator_t>                          MyMovableShmMultiSet;
typedef map<test::movable_int, test::movable_int, 
      std::less<test::movable_int>, 
      shmem_movable_node_pair_allocator_t>                     MyMovableShmMap;
typedef multimap<test::movable_int, test::movable_int, 
      std::less<test::movable_int>, 
      shmem_movable_node_pair_allocator_t>                     MyMovableShmMultiMap;

typedef set<test::movable_and_copyable_int
           ,std::less<test::movable_and_copyable_int>
           ,shmem_move_copy_allocator_t>                  MyMoveCopyShmSet;
typedef multiset<test::movable_and_copyable_int, 
      std::less<test::movable_and_copyable_int>, 
      shmem_move_copy_allocator_t>                        MyMoveCopyShmMultiSet;
typedef map<test::movable_and_copyable_int
           ,test::movable_and_copyable_int
           ,std::less<test::movable_and_copyable_int>
           ,shmem_move_copy_node_pair_allocator_t>             MyMoveCopyShmMap;
typedef multimap<test::movable_and_copyable_int
                ,test::movable_and_copyable_int
                ,std::less<test::movable_and_copyable_int>
                ,shmem_move_copy_node_pair_allocator_t>        MyMoveCopyShmMultiMap;
//Test recursive structures
class recursive_set
{
public:
   int id_;
   set<recursive_set> set_;
   friend bool operator< (const recursive_set &a, const recursive_set &b)
   {  return a.id_ < b.id_;   }
};

class recursive_map
{
   public:
   int id_;
   map<recursive_map, recursive_map> map_;
   friend bool operator< (const recursive_map &a, const recursive_map &b)
   {  return a.id_ < b.id_;   }
};

//Test recursive structures
class recursive_multiset
{
public:
   int id_;
   multiset<recursive_multiset> multiset_;
   friend bool operator< (const recursive_multiset &a, const recursive_multiset &b)
   {  return a.id_ < b.id_;   }
};

class recursive_multimap
{
public:
   int id_;
   multimap<recursive_multimap, recursive_multimap> multimap_;
   friend bool operator< (const recursive_multimap &a, const recursive_multimap &b)
   {  return a.id_ < b.id_;   }
};

template<class C>
void test_move()
{
   //Now test move semantics
   C original;
   C move_ctor(boost::interprocess::move(original));
   C move_assign;
   move_assign = boost::interprocess::move(move_ctor);
   move_assign.swap(original);
}

int main ()
{
   //Recursive container instantiation
   {
      set<recursive_set> set_;
      multiset<recursive_multiset> multiset_;
      map<recursive_map, recursive_map> map_;
      multimap<recursive_multimap, recursive_multimap> multimap_;
   }
   //Now test move semantics
   {
      test_move<set<recursive_set> >();
      test_move<multiset<recursive_multiset> >();
      test_move<map<recursive_map, recursive_map> >();
      test_move<multimap<recursive_multimap, recursive_multimap> >();
   }

   using namespace boost::interprocess::detail;

   if(0 != test::set_test<my_managed_shared_memory
                        ,MyShmSet
                        ,MyStdSet
                        ,MyShmMultiSet
                        ,MyStdMultiSet>()){
      return 1;
   }

   if(0 != test::set_test_copyable<my_managed_shared_memory
                        ,MyShmSet
                        ,MyStdSet
                        ,MyShmMultiSet
                        ,MyStdMultiSet>()){
      return 1;
   }

   if(0 != test::set_test<my_managed_shared_memory
                        ,MyMovableShmSet
                        ,MyStdSet
                        ,MyMovableShmMultiSet
                        ,MyStdMultiSet>()){
      return 1;
   }

   if(0 != test::set_test<my_managed_shared_memory
                        ,MyMoveCopyShmSet
                        ,MyStdSet
                        ,MyMoveCopyShmMultiSet
                        ,MyStdMultiSet>()){
      return 1;
   }


   if (0 != test::map_test<my_managed_shared_memory
                  ,MyShmMap
                  ,MyStdMap
                  ,MyShmMultiMap
                  ,MyStdMultiMap>()){
      return 1;
   }

   if(0 != test::map_test_copyable<my_managed_shared_memory
                        ,MyShmMap
                        ,MyStdMap
                        ,MyShmMultiMap
                        ,MyStdMultiMap>()){
      return 1;
   }

//   if (0 != test::map_test<my_managed_shared_memory
//                  ,MyMovableShmMap
//                  ,MyStdMap
//                  ,MyMovableShmMultiMap
//                  ,MyStdMultiMap>()){
//      return 1;
//   }

   if (0 != test::map_test<my_managed_shared_memory
                  ,MyMoveCopyShmMap
                  ,MyStdMap
                  ,MyMoveCopyShmMultiMap
                  ,MyStdMultiMap>()){
      return 1;
   }

   const test::EmplaceOptions SetOptions = (test::EmplaceOptions)(test::EMPLACE_HINT | test::EMPLACE_ASSOC);
   if(!boost::interprocess::test::test_emplace<set<test::EmplaceInt>, SetOptions>())
      return 1;
   if(!boost::interprocess::test::test_emplace<multiset<test::EmplaceInt>, SetOptions>())
      return 1;
   const test::EmplaceOptions MapOptions = (test::EmplaceOptions)(test::EMPLACE_HINT_PAIR | test::EMPLACE_ASSOC_PAIR);
   if(!boost::interprocess::test::test_emplace<map<test::EmplaceInt, test::EmplaceInt>, MapOptions>())
      return 1;
   if(!boost::interprocess::test::test_emplace<multimap<test::EmplaceInt, test::EmplaceInt>, MapOptions>())
      return 1;

   return 0;
}

#include <boost/interprocess/detail/config_end.hpp>