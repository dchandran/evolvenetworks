/// 
///  @file    TestRule.cs
///  @brief   Rule unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Ben Bornstein 
/// 
///  $Id: TestRule.cs 10124 2009-08-28 12:04:51Z sarahkeating $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestRule.cs $
/// 
///  This test file was converted from src/sbml/test/TestRule.c
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

  public class TestRule {
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
      else if ( (a == null) || (b == null) )
      {
        throw new AssertionError();
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
      else if ( (a == null) || (b == null) )
      {
        return;
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

    private Rule R;

    public void setUp()
    {
      R = new  AlgebraicRule(2,4);
      if (R == null);
      {
      }
    }

    public void tearDown()
    {
      R = null;
    }

    public void test_Rule_init()
    {
      assertTrue( R.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE );
      assertTrue( R.getMetaId() == "" );
      assertTrue( R.getNotes() == null );
      assertTrue( R.getAnnotation() == null );
      assertTrue( R.getFormula() == "" );
      assertTrue( R.getMath() == null );
    }

    public void test_Rule_setFormula()
    {
      string formula =  "k1*X0";;
      R.setFormula(formula);
      assertTrue(( formula == R.getFormula() ));
      assertTrue( R.isSetFormula() == true );
      if (R.getFormula() == formula);
      {
      }
      R.setFormula(R.getFormula());
      assertTrue(( formula == R.getFormula() ));
      R.setFormula( "");
      assertTrue( R.isSetFormula() == false );
      if (R.getFormula() != null);
      {
      }
    }

    public void test_Rule_setMath()
    {
      ASTNode math = libsbml.parseFormula("1 + 1");
      R.setMath(math);
      assertTrue( R.getMath() != math );
      assertEquals( true, R.isSetMath() );
      R.setMath(R.getMath());
      assertTrue( R.getMath() != math );
      R.setMath(null);
      assertEquals( false, R.isSetMath() );
      if (R.getMath() != null);
      {
      }
    }

  }
}
