#
# @file    TestReadFromFile8.rb
# @brief   Reads test-data/l2v4-new.xml into memory and tests it.
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Sarah Keating 
#
# $Id$
# $HeadURL$
#
# This test file was converted from src/sbml/test/TestReadFromFile8.cpp
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
require 'test/unit'
require 'libSBML'

class TestReadFromFile8 < Test::Unit::TestCase

  def test_read_l2v4_new
    reader = LibSBML::SBMLReader.new()
    filename = "../../sbml/test/test-data/"
    filename += "l2v4-new.xml"
    d = reader.readSBML(filename)
    if (d == nil)
    end
    assert( d.getLevel() == 2 )
    assert( d.getVersion() == 4 )
    m = d.getModel()
    assert( m != nil )
    assert( m.getId() ==  "l2v4_all" )
    assert( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    assert( c != nil )
    assert( c.getId() ==  "a" )
    assert( c.getSize() == 1 )
    assert_equal false, c.getConstant()
    assert( m.getNumEvents() == 1 )
    e = m.getEvent(0)
    assert( e != nil )
    assert_equal true, e.getUseValuesFromTriggerTime()
    assert_equal true, e.isSetTrigger()
    trigger = e.getTrigger()
    assert( trigger != nil )
    ast = trigger.getMath()
    assert ((  "lt(x, 3)" == LibSBML::formulaToString(ast) ))
    assert( e.getNumEventAssignments() == 1 )
    ea = e.getEventAssignment(0)
    assert( ea != nil )
    assert( ea.getVariable() ==  "a" )
    ast = ea.getMath()
    assert ((  "x * p3" == LibSBML::formulaToString(ast) ))
    d = nil
  end

end
