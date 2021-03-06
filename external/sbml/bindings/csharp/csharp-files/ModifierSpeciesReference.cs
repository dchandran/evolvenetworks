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

public class ModifierSpeciesReference : SimpleSpeciesReference {
	private HandleRef swigCPtr;
	
	internal ModifierSpeciesReference(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.ModifierSpeciesReferenceUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.ModifierSpeciesReferenceUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(ModifierSpeciesReference obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (ModifierSpeciesReference obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~ModifierSpeciesReference() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_ModifierSpeciesReference(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public ModifierSpeciesReference(long level, long version) : this(libsbmlPINVOKE.new_ModifierSpeciesReference__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public ModifierSpeciesReference(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_ModifierSpeciesReference__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public override SBase clone() {
    IntPtr cPtr = libsbmlPINVOKE.ModifierSpeciesReference_clone(swigCPtr);
    ModifierSpeciesReference ret = (cPtr == IntPtr.Zero) ? null : new ModifierSpeciesReference(cPtr, true);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.ModifierSpeciesReference_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.ModifierSpeciesReference_getElementName(swigCPtr);
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.ModifierSpeciesReference_hasRequiredAttributes(swigCPtr);
    return ret;
  }

}

}
