#include "RoadkillCmd.h"

#include <maya/MFnPlugin.h>
#include <maya/MIOStream.h>
#include <maya/MGlobal.h>

#include <time.h>


#ifdef WIN32
#define EXTERN_DECL __declspec( dllexport )
#else
#define EXTERN_DECL extern
#endif


char	gErrorMessage[512] = {0};


EXTERN_DECL MStatus initializePlugin( MObject obj );
EXTERN_DECL MStatus uninitializePlugin( MObject obj );



/****************************************
*		Initialise the plug-in			*
****************************************/

MStatus initializePlugin(MObject obj)
{
	MStatus status;

	MFnPlugin plugin( obj, "(c) 2009 - 2018 Andy Swann", "Version 1.045 R+D", "Any");

	status = plugin.registerCommand("RoadkillPro", RoadkillCmd::creator, RoadkillCmd::newSyntax);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MGlobal::executeCommand("RoadkillProMenu", false, false);
	
	return MS::kSuccess;
}



/****************************************
*	   Uninitialise the plug-in			*
****************************************/

MStatus uninitializePlugin(MObject obj)
{
	MStatus status;

	MFnPlugin plugin( obj );

	status = plugin.deregisterCommand("RoadkillPro");
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return MS::kSuccess;
}
