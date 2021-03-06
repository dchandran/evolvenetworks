#
# @file    TestCompartment_newSetters.py
# @brief   Compartment unit tests for new set function API
#
# @author  Akiya Jouraku (Python conversion)
# @author  Sarah Keating 
#
# $Id$
# $HeadURL$
#
# This test file was converted from src/sbml/test/TestCompartment_newSetters.c
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

class TestCompartment_newSetters(unittest.TestCase):

  C = None

  def setUp(self):
    self.C = libsbml.Compartment(1,2)
    if (self.C == None):
      pass    
    pass  

  def tearDown(self):
    self.C = None
    pass  

  def test_Compartment_setCompartmentType1(self):
    i = self.C.setCompartmentType( "cell")
    self.assert_( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE )
    self.assertEqual( False, self.C.isSetCompartmentType() )
    i = self.C.unsetCompartmentType()
    self.assert_( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE )
    self.assertEqual( False, self.C.isSetCompartmentType() )
    pass  

  def test_Compartment_setCompartmentType2(self):
    c = libsbml.Compartment(2,2)
    i = c.setCompartmentType( "1cell")
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assertEqual( False, c.isSetCompartmentType() )
    i = c.unsetCompartmentType()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetCompartmentType() )
    c = None
    pass  

  def test_Compartment_setCompartmentType3(self):
    c = libsbml.Compartment(2,2)
    i = c.setCompartmentType( "cell")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( True, c.isSetCompartmentType() )
    self.assert_((  "cell"  == c.getCompartmentType() ))
    i = c.unsetCompartmentType()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetCompartmentType() )
    c = None
    pass  

  def test_Compartment_setCompartmentType4(self):
    c = libsbml.Compartment(2,2)
    i = c.setCompartmentType("")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetCompartmentType() )
    c = None
    pass  

  def test_Compartment_setConstant1(self):
    i = self.C.setConstant(False)
    self.assert_( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE )
    self.assert_( self.C.getConstant() == False )
    pass  

  def test_Compartment_setConstant2(self):
    c = libsbml.Compartment(2,2)
    i = c.setConstant(False)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( c.getConstant() == False )
    c = None
    pass  

  def test_Compartment_setId2(self):
    c = libsbml.Compartment(2,2)
    i = c.setId( "1cell")
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assertEqual( False, c.isSetId() )
    c = None
    pass  

  def test_Compartment_setId3(self):
    c = libsbml.Compartment(2,2)
    i = c.setId( "cell")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( True, c.isSetId() )
    self.assert_((  "cell"  == c.getId() ))
    c = None
    pass  

  def test_Compartment_setId4(self):
    c = libsbml.Compartment(2,2)
    i = c.setId("")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetId() )
    c = None
    pass  

  def test_Compartment_setName1(self):
    i = self.C.setName( "cell")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( True, self.C.isSetName() )
    i = self.C.unsetName()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetName() )
    pass  

  def test_Compartment_setName2(self):
    c = libsbml.Compartment(1,2)
    i = c.setName( "1cell")
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assertEqual( False, c.isSetName() )
    i = c.unsetName()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetName() )
    c = None
    pass  

  def test_Compartment_setName3(self):
    i = self.C.setName("")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetName() )
    pass  

  def test_Compartment_setOutside1(self):
    i = self.C.setOutside( "1cell")
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assertEqual( False, self.C.isSetOutside() )
    i = self.C.unsetOutside()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetOutside() )
    pass  

  def test_Compartment_setOutside2(self):
    i = self.C.setOutside( "litre")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( True, self.C.isSetOutside() )
    i = self.C.unsetOutside()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetOutside() )
    pass  

  def test_Compartment_setOutside3(self):
    i = self.C.setOutside("")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetOutside() )
    pass  

  def test_Compartment_setSize1(self):
    i = self.C.setSize(2.0)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( self.C.getSize() == 2.0 )
    i = self.C.unsetSize()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    pass  

  def test_Compartment_setSize2(self):
    c = libsbml.Compartment(2,2)
    i = c.setSize(4)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( c.getSize() == 4 )
    self.assertEqual( True, c.isSetSize() )
    i = c.unsetSize()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetSize() )
    c = None
    pass  

  def test_Compartment_setSpatialDimensions1(self):
    i = self.C.setSpatialDimensions(2)
    self.assert_( i == libsbml.LIBSBML_UNEXPECTED_ATTRIBUTE )
    self.assert_( self.C.getSpatialDimensions() == 3 )
    pass  

  def test_Compartment_setSpatialDimensions2(self):
    c = libsbml.Compartment(2,2)
    i = c.setSpatialDimensions(4)
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assert_( c.getSpatialDimensions() == 3 )
    c = None
    pass  

  def test_Compartment_setSpatialDimensions3(self):
    c = libsbml.Compartment(2,2)
    i = c.setSpatialDimensions(2)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( c.getSpatialDimensions() == 2 )
    c = None
    pass  

  def test_Compartment_setUnits1(self):
    i = self.C.setUnits( "1cell")
    self.assert_( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE )
    self.assertEqual( False, self.C.isSetUnits() )
    i = self.C.unsetUnits()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetUnits() )
    pass  

  def test_Compartment_setUnits2(self):
    i = self.C.setUnits( "litre")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( True, self.C.isSetUnits() )
    i = self.C.unsetUnits()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetUnits() )
    pass  

  def test_Compartment_setUnits3(self):
    i = self.C.setUnits("")
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, self.C.isSetUnits() )
    pass  

  def test_Compartment_setVolume1(self):
    i = self.C.setVolume(2.0)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( self.C.getVolume() == 2.0 )
    self.assertEqual( True, self.C.isSetVolume() )
    i = self.C.unsetVolume()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( self.C.getVolume() == 1.0 )
    self.assertEqual( True, self.C.isSetVolume() )
    pass  

  def test_Compartment_setVolume2(self):
    c = libsbml.Compartment(2,2)
    i = c.setVolume(4)
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assert_( c.getVolume() == 4.0 )
    self.assertEqual( True, c.isSetVolume() )
    i = c.unsetVolume()
    self.assert_( i == libsbml.LIBSBML_OPERATION_SUCCESS )
    self.assertEqual( False, c.isSetVolume() )
    c = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestCompartment_newSetters))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
