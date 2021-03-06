#
# @file    TestConsistencyChecks.py
# @brief   Reads test-data/inconsistent.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Sarah Keating 
#
# $Id: TestConsistencyChecks.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/sbml/TestConsistencyChecks.py $
#
# This test file was converted from src/sbml/test/TestConsistencyChecks.cpp
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

class TestConsistencyChecks(unittest.TestCase):


  def test_consistency_checks(self):
    reader = libsbml.SBMLReader()
    d = libsbml.SBMLDocument()
    filename = "../../sbml/test/test-data/"
    filename += "inconsistent.xml"
    d = reader.readSBML(filename)
    if (d == None):
      pass    
    errors = d.checkConsistency()
    self.assert_( errors == 1 )
    self.assert_( d.getError(0).getErrorId() == 10301 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY,False)
    errors = d.checkConsistency()
    self.assert_( errors == 1 )
    self.assert_( d.getError(0).getErrorId() == 20612 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY,False)
    errors = d.checkConsistency()
    self.assert_( errors == 1 )
    self.assert_( d.getError(0).getErrorId() == 10701 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY,False)
    errors = d.checkConsistency()
    self.assert_( errors == 1 )
    self.assert_( d.getError(0).getErrorId() == 10214 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY,False)
    errors = d.checkConsistency()
    self.assert_( errors == 2 )
    self.assert_( d.getError(0).getErrorId() == 10523 )
    self.assert_( d.getError(1).getErrorId() == 99505 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,False)
    errors = d.checkConsistency()
    self.assert_( errors == 0 )
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestConsistencyChecks))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
