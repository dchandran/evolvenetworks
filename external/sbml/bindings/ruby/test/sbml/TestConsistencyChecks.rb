#
# @file    TestConsistencyChecks.rb
# @brief   Reads test-data/inconsistent.xml into memory and tests it.
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Sarah Keating 
#
# $Id: TestConsistencyChecks.rb 10124 2009-08-28 12:04:51Z sarahkeating $
# $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/ruby/test/sbml/TestConsistencyChecks.rb $
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
require 'test/unit'
require 'libSBML'

class TestConsistencyChecks < Test::Unit::TestCase

  def test_consistency_checks
    reader = LibSBML::SBMLReader.new()
    d = LibSBML::SBMLDocument.new()
    filename = "../../sbml/test/test-data/"
    filename += "inconsistent.xml"
    d = reader.readSBML(filename)
    if (d == nil)
    end
    errors = d.checkConsistency()
    assert( errors == 1 )
    assert( d.getError(0).getErrorId() == 10301 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(LibSBML::LIBSBML_CAT_IDENTIFIER_CONSISTENCY,false)
    errors = d.checkConsistency()
    assert( errors == 1 )
    assert( d.getError(0).getErrorId() == 20612 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(LibSBML::LIBSBML_CAT_GENERAL_CONSISTENCY,false)
    errors = d.checkConsistency()
    assert( errors == 1 )
    assert( d.getError(0).getErrorId() == 10701 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(LibSBML::LIBSBML_CAT_SBO_CONSISTENCY,false)
    errors = d.checkConsistency()
    assert( errors == 1 )
    assert( d.getError(0).getErrorId() == 10214 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(LibSBML::LIBSBML_CAT_MATHML_CONSISTENCY,false)
    errors = d.checkConsistency()
    assert( errors == 2 )
    assert( d.getError(0).getErrorId() == 10523 )
    assert( d.getError(1).getErrorId() == 99505 )
    d.getErrorLog().clearLog()
    d.setConsistencyChecks(LibSBML::LIBSBML_CAT_UNITS_CONSISTENCY,false)
    errors = d.checkConsistency()
    assert( errors == 0 )
    d = nil
  end

end
