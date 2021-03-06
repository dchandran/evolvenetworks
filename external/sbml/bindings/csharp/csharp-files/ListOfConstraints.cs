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

public class ListOfConstraints : ListOf {
	private HandleRef swigCPtr;
	
	internal ListOfConstraints(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.ListOfConstraintsUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.ListOfConstraintsUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(ListOfConstraints obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (ListOfConstraints obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~ListOfConstraints() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_ListOfConstraints(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public new ListOfConstraints clone() {
    IntPtr cPtr = libsbmlPINVOKE.ListOfConstraints_clone(swigCPtr);
    ListOfConstraints ret = (cPtr == IntPtr.Zero) ? null : new ListOfConstraints(cPtr, true);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.ListOfConstraints_getTypeCode(swigCPtr);
    return ret;
  }

  public override int getItemTypeCode() {
    int ret = libsbmlPINVOKE.ListOfConstraints_getItemTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.ListOfConstraints_getElementName(swigCPtr);
    return ret;
  }

  public new Constraint get(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfConstraints_get__SWIG_0(swigCPtr, n);
    Constraint ret = (cPtr == IntPtr.Zero) ? null : new Constraint(cPtr, false);
    return ret;
  }

  public new Constraint remove(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfConstraints_remove(swigCPtr, n);
    Constraint ret = (cPtr == IntPtr.Zero) ? null : new Constraint(cPtr, true);
    return ret;
  }

  public ListOfConstraints() : this(libsbmlPINVOKE.new_ListOfConstraints(), true) {
  }

}

}
