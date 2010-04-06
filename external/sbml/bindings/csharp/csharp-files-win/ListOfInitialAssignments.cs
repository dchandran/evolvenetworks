/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

namespace libsbml {

using System;
using System.Runtime.InteropServices;

public class ListOfInitialAssignments : ListOf {
	private HandleRef swigCPtr;
	
	internal ListOfInitialAssignments(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.ListOfInitialAssignmentsUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.ListOfInitialAssignmentsUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(ListOfInitialAssignments obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (ListOfInitialAssignments obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~ListOfInitialAssignments() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_ListOfInitialAssignments(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public new ListOfInitialAssignments clone() {
    IntPtr cPtr = libsbmlPINVOKE.ListOfInitialAssignments_clone(swigCPtr);
    ListOfInitialAssignments ret = (cPtr == IntPtr.Zero) ? null : new ListOfInitialAssignments(cPtr, true);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.ListOfInitialAssignments_getTypeCode(swigCPtr);
    return ret;
  }

  public override int getItemTypeCode() {
    int ret = libsbmlPINVOKE.ListOfInitialAssignments_getItemTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.ListOfInitialAssignments_getElementName(swigCPtr);
    return ret;
  }

  public new InitialAssignment get(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfInitialAssignments_get__SWIG_0(swigCPtr, n);
    InitialAssignment ret = (cPtr == IntPtr.Zero) ? null : new InitialAssignment(cPtr, false);
    return ret;
  }

  public new InitialAssignment get(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfInitialAssignments_get__SWIG_2(swigCPtr, sid);
    InitialAssignment ret = (cPtr == IntPtr.Zero) ? null : new InitialAssignment(cPtr, false);
    return ret;
  }

  public new InitialAssignment remove(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfInitialAssignments_remove__SWIG_0(swigCPtr, n);
    InitialAssignment ret = (cPtr == IntPtr.Zero) ? null : new InitialAssignment(cPtr, true);
    return ret;
  }

  public new InitialAssignment remove(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfInitialAssignments_remove__SWIG_1(swigCPtr, sid);
    InitialAssignment ret = (cPtr == IntPtr.Zero) ? null : new InitialAssignment(cPtr, true);
    return ret;
  }

  public ListOfInitialAssignments() : this(libsbmlPINVOKE.new_ListOfInitialAssignments(), true) {
  }

}

}
