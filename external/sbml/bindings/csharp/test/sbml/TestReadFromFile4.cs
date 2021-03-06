/// 
///  @file    TestReadFromFile4.cs
///  @brief   Reads tests/l1v1-minimal.xml into memory and tests it.
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Ben Bornstein 
/// 
///  $Id: TestReadFromFile4.cs 10124 2009-08-28 12:04:51Z sarahkeating $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestReadFromFile4.cs $
/// 
///  This test file was converted from src/sbml/test/TestReadFromFile4.c
///  with the help of conversion sciprt (ctest_converter.pl).
/// 
/// <!---------------------------------------------------------------------------
///  This file is part of libSBML.  Please visit http://sbml.org for more
///  information about SBML, and the latest version of libSBML.
/// 
///  Copyright 2005-2009 California Institute of Technology.
///  Copyright 2002-2005 California Institute of Technology and
///                      Japan Science and Technology Corporation.
///  
///  This library is free software; you can redistribute it and/or modify it
///  under the terms of the GNU Lesser General Public License as published by
///  the Free Software Foundation.  A copy of the license agreement is provided
///  in the file named "LICENSE.txt" included with this software distribution
///  and also available online as http://sbml.org/software/libsbml/license.html
/// --------------------------------------------------------------------------->*/


namespace LibSBMLCSTest {

  using libsbml;

  using  System.IO;

  public class TestReadFromFile4 {
    public class AssertionError : System.Exception 
    {
      public AssertionError() : base()
      {
        
      }
    }


    static void assertTrue(bool condition)
    {
      if (condition == true)
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        return;
      }
      else if (a.Equals(b))
      {
        return;
      }
  
      throw new AssertionError();
    }

    static void assertNotEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        throw new AssertionError();
      }
      else if (a.Equals(b))
      {
        throw new AssertionError();
      }
    }

    static void assertEquals(bool a, bool b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(bool a, bool b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(int a, int b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(int a, int b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }


    public void test_read_l1v1_minimal()
    {
      SBMLDocument d;
      Model m;
      Compartment c;
      Reaction r;
      Species s;
      SpeciesReference sr;
      string filename = "../../sbml/test/test-data/l1v1-minimal.xml";
      d = libsbml.readSBML(filename);
      if (d == null);
      {
      }
      assertTrue( d.getLevel() == 1 );
      assertTrue( d.getVersion() == 1 );
      m = d.getModel();
      assertTrue( m.getNumCompartments() == 1 );
      c = m.getCompartment(0);
      assertTrue((  "x" == c.getName() ));
      assertTrue( m.getNumSpecies() == 1 );
      s = m.getSpecies(0);
      assertTrue((  "y"  == s.getName() ));
      assertTrue((  "x"  == s.getCompartment() ));
      assertTrue( s.getInitialAmount() == 1 );
      assertTrue( s.getBoundaryCondition() == false );
      assertTrue( m.getNumReactions() == 1 );
      r = m.getReaction(0);
      assertTrue((  "r" == r.getName() ));
      assertTrue( r.getReversible() != false );
      assertTrue( r.getFast() == false );
      assertTrue( r.getNumReactants() == 1 );
      assertTrue( r.getNumProducts() == 1 );
      sr = r.getReactant(0);
      assertTrue((  "y" == sr.getSpecies() ));
      assertTrue( sr.getStoichiometry() == 1 );
      assertTrue( sr.getDenominator() == 1 );
      sr = r.getProduct(0);
      assertTrue((  "y" == sr.getSpecies() ));
      assertTrue( sr.getStoichiometry() == 1 );
      assertTrue( sr.getDenominator() == 1 );
      d = null;
    }

  }
}
