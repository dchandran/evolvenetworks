#
# @file    TestListOf.py
# @brief   ListOf unit tests
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id: TestListOf.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestListOf.py $
#
# This test file was converted from src/sbml/test/TestListOf.c
# with the help of conversion sciprt (ctest_converter.pl).
#
#<!---------------------------------------------------------------------------
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright 2005-2009 California Institute of Technology.
# Copyright 2002-2005 California Institute of Technology and
#                     Japan Science and Technology Corporation.
# 
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
#--------------------------------------------------------------------------->*/
import sys
import unittest
import libsbml

class TestListOf(unittest.TestCase):

  def test_ListOf_clear(self):
    lo = libsbml.ListOf()
    sp = libsbml.Species(2,4)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    self.assert_( lo.size() == 5 )
    lo.clear(1)
    self.assert_( lo.size() == 0 )
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.appendAndOwn(sp)
    self.assert_( lo.size() == 5 )
    lo.get(0)
    lo.get(1)
    lo.get(2)
    lo.get(3)
    lo.get(4)
    lo.clear(0)
    self.assert_( lo.size() == 0 )
    lo = None
    pass  

  def test_ListOf_create(self):
    lo = libsbml.ListOf()
    self.assert_( lo.getTypeCode() == libsbml.SBML_LIST_OF )
    self.assert_( lo.getNotes() == None )
    self.assert_( lo.getAnnotation() == None )
    self.assert_( lo.getMetaId() == "" )
    self.assert_( lo.size() == 0 )
    lo = None
    pass  

  def test_ListOf_free_NULL(self):
    pass  

  def test_ListOf_remove(self):
    lo = libsbml.ListOf()
    sp = libsbml.Species(2,4)
    self.assert_( lo.size() == 0 )
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    self.assert_( lo.size() == 5 )
    lo.remove(0)
    lo.remove(0)
    lo.remove(0)
    lo.remove(0)
    lo.remove(0)
    self.assert_( lo.size() == 0 )
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.append(sp)
    lo.appendAndOwn(sp)
    self.assert_( lo.size() == 5 )
    lo = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestListOf))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
