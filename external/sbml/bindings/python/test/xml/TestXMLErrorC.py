#
# @file    TestXMLErrorC.py
# @brief   XMLError unit tests, C version
#
# @author  Akiya Jouraku (Python conversion)
# @author  Sarah Keating 
#
# $Id: TestXMLErrorC.py 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test/xml/TestXMLErrorC.py $
#
# This test file was converted from src/sbml/test/TestXMLErrorC.c
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

class TestXMLErrorC(unittest.TestCase):


  def test_XMLError_create_C(self):
    error = libsbml.XMLError()
    self.assert_( error != None )
    self.assert_( error.isInfo() == False )
    self.assert_( error.isWarning() == False )
    self.assert_( error.isError() == False )
    self.assert_( error.isFatal() == True )
    error = None
    error = libsbml.XMLError(12345, "My message")
    self.assert_( (  "My message" != error.getMessage() ) == False )
    self.assert_( error.getErrorId() == 12345 )
    error = None
    pass  

  def test_XMLError_variablesAsStrings(self):
    error = libsbml.XMLError(1003, "")
    self.assert_( error.getErrorId() == 1003 )
    self.assert_( error.getSeverity() == libsbml.LIBSBML_SEV_ERROR )
    self.assert_((  "Error" == error.getSeverityAsString() ))
    self.assert_( error.getCategory() == libsbml.LIBSBML_CAT_XML )
    self.assert_((  "XML content" == error.getCategoryAsString() ))
    error = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestXMLErrorC))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
