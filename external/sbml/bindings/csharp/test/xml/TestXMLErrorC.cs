/// 
///  @file    TestXMLErrorC.cs
///  @brief   XMLError unit tests, C version
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestXMLErrorC.cs 10124 2009-08-28 12:04:51Z sarahkeating $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/xml/TestXMLErrorC.cs $
/// 
///  This test file was converted from src/sbml/test/TestXMLErrorC.c
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

  public class TestXMLErrorC {
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


    public void test_XMLError_create_C()
    {
      XMLError error = new  XMLError();
      assertTrue( error != null );
      assertTrue( error.isInfo() == false );
      assertTrue( error.isWarning() == false );
      assertTrue( error.isError() == false );
      assertTrue( error.isFatal() == true );
      error = null;
      error = new  XMLError(12345, "My message");
      assertTrue( (  "My message" != error.getMessage() ) == false );
      assertTrue( error.getErrorId() == 12345 );
      error = null;
    }

    public void test_XMLError_variablesAsStrings()
    {
      XMLError error = new  XMLError(1003, "");
      assertTrue( error.getErrorId() == 1003 );
      assertTrue( error.getSeverity() == libsbml.LIBSBML_SEV_ERROR );
      assertTrue((  "Error" == error.getSeverityAsString() ));
      assertTrue( error.getCategory() == libsbml.LIBSBML_CAT_XML );
      assertTrue((  "XML content" == error.getCategoryAsString() ));
      error = null;
    }

  }
}
