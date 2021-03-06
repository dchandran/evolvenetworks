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

public class SpeciesType : SBase {
	private HandleRef swigCPtr;
	
	internal SpeciesType(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.SpeciesTypeUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.SpeciesTypeUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(SpeciesType obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (SpeciesType obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~SpeciesType() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_SpeciesType(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public SpeciesType(long level, long version) : this(libsbmlPINVOKE.new_SpeciesType__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public SpeciesType(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_SpeciesType__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public SpeciesType(SpeciesType orig) : this(libsbmlPINVOKE.new_SpeciesType__SWIG_2(SpeciesType.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public new SpeciesType clone() {
    IntPtr cPtr = libsbmlPINVOKE.SpeciesType_clone(swigCPtr);
    SpeciesType ret = (cPtr == IntPtr.Zero) ? null : new SpeciesType(cPtr, true);
    return ret;
  }

  public new string getId() {
    string ret = libsbmlPINVOKE.SpeciesType_getId(swigCPtr);
    return ret;
  }

  public new string getName() {
    string ret = libsbmlPINVOKE.SpeciesType_getName(swigCPtr);
    return ret;
  }

  public new bool isSetId() {
    bool ret = libsbmlPINVOKE.SpeciesType_isSetId(swigCPtr);
    return ret;
  }

  public new bool isSetName() {
    bool ret = libsbmlPINVOKE.SpeciesType_isSetName(swigCPtr);
    return ret;
  }

  public new int setId(string sid) {
    int ret = libsbmlPINVOKE.SpeciesType_setId(swigCPtr, sid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public new int setName(string name) {
    int ret = libsbmlPINVOKE.SpeciesType_setName(swigCPtr, name);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int unsetName() {
    int ret = libsbmlPINVOKE.SpeciesType_unsetName(swigCPtr);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.SpeciesType_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.SpeciesType_getElementName(swigCPtr);
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.SpeciesType_hasRequiredAttributes(swigCPtr);
    return ret;
  }

}

}
