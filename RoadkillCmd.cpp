#include "RoadkillCmd.h"

#include "Unwrapper/RKProgress.h"
#include "Unwrapper/Unwrap.h"

#include <algorithm>
#include <math.h>

#include <maya/MSyntax.h>
#include <maya/MDagPath.h>
#include <maya/MIOStream.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnTransform.h>

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MFnSingleIndexedComponent.h>




extern char gErrorMessage[512];

extern int gFunction;
extern int gScale;
extern int gUpright;
bool gAllInOnePack = false;

bool gSetupFailed = false;

extern int gMaintainUVs;
extern int gMatch;
extern int gEdgeLen;
extern int gSymMatch;
extern double gMatchAngleTolerence;
extern double gMatchEdgeTolerence;

extern int gOuterBorderSwitch;
extern int gOuterBorderPixels;
extern int gOuterBorderTMapSize;





/****************************************
*			Command Constructor			*
****************************************/

RoadkillCmd::RoadkillCmd()
{
	m_WorkingList.clear();
}


/****************************************
*		  Command Deconstructor			*
****************************************/

RoadkillCmd::~RoadkillCmd() 
{
	ClearWorkingList();
}


bool RoadkillCmd::isUndoable() const
{
	if(gSetupFailed == true) return(false);

	return(true);
}



/****************************************
*		  Command Creator				*
*	Split all dodgy edges here for		*
*	Every selected mesh, then reset		*
*	The selection list					*
****************************************/

void* RoadkillCmd::creator() 
{ 
	return new RoadkillCmd;
}




void RoadkillCmd::Setup()
{
	gSetupFailed = false;


	MSelectionList selList;
	MGlobal::getActiveSelectionList(selList);
	MItSelectionList selListIter(selList);
	selListIter.setFilter(MFn::kMesh);

	RKProgress::Get().NewProgressWindow();
	RKProgress::Get().SetText("Initialising");
	RKProgress::Get().SetProgress(0);

	bool foundShite = false;
	for( ; !selListIter.isDone(); selListIter.next() )
	{
		MDagPath dagPath;
		MObject component;
		selListIter.getDagPath(dagPath, component);

		if(dagPath.apiType() == MFn::kMesh)
		{
			WorkingMesh* pNewMesh = new WorkingMesh();
			pNewMesh->m_MeshNode = dagPath;

			MFn::Type testType = component.apiType();
			if(testType != MFn::kInvalid)
			{
				foundShite = true;					// found something we don't want to deal with

				if(component.apiType() == MFn::kMeshEdgeComponent)
				{
					MFnSingleIndexedComponent compFn(component);
					compFn.getElements(pNewMesh->m_SelectedEdges);
					foundShite = false;
				}

				if(component.apiType() == MFn::kMeshPolygonComponent)
				{
					MFnSingleIndexedComponent compFn(component);
					compFn.getElements(pNewMesh->m_SelectedFaces);
					foundShite = false;
				}
			}

			m_WorkingList.push_back(pNewMesh);
		}
	}


	if(m_WorkingList.size() == 0 || foundShite)
	{
		displayError( "Roadkill: Please select mesh(es) or mesh faces" );
		RKProgress::Get().CloseProgressWindow();
		gSetupFailed = true;
		return;
	}


	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		m_WorkingList[Index]->SplitSingletons();

		if(m_WorkingList[Index]->m_untextured == true)
		{
			displayError( "Roadkill: Attempting to work with a mesh or faces with no UVs" );
			RKProgress::Get().CloseProgressWindow();
			gSetupFailed = true;
			return;
		}
	}


	MSelectionList StoreHiliteList;
	MSelectionList StoreSelList;


	if(gFunction == 6 || gFunction == 8)
	{
		if(gFunction == 6)
		{
			if(m_WorkingList.size() < 2)
			{
				displayError( "Roadkill: Only one mesh selected." );
			}

			ScaleToSource();
		}
		else
		{
			SymmetryCut();
		}
	}
	else
	{
		// don't bother with this for packing, or scaling

		MGlobal::getActiveSelectionList(StoreSelList);

		for(int Index = 0; Index < m_WorkingList.size(); Index++)
		{
			m_WorkingList[Index]->issuePolyCut(true);				// find history!
		}

		MGlobal::setActiveSelectionList(StoreSelList);

		Unwrap* pUnwrapper;

		if(gAllInOnePack)
		{
			pUnwrapper = new Unwrap();

			for(int Index = 0; Index < m_WorkingList.size(); Index++)
			{
				m_WorkingList[Index]->AddToUnwrapper(pUnwrapper);
			}

			RKProgress::Get().SetText("Unwrapping");
			RKProgress::Get().SetProgress(0);
			pUnwrapper->DoUnwrap(true);
			RKProgress::Get().SetProgress(100);

			for(int Index = 0; Index < m_WorkingList.size(); Index++)
			{
				m_WorkingList[Index]->SetNewUVs(pUnwrapper);
			}

			pUnwrapper->DeleteVertsAndFaces();
			delete pUnwrapper;
		}
		else
		{
			for(int Index = 0; Index < m_WorkingList.size(); Index++)
			{
				pUnwrapper = new Unwrap();
				m_WorkingList[Index]->AddToUnwrapper(pUnwrapper);

				RKProgress::Get().SetText("Unwrapping");
				RKProgress::Get().SetProgress(0);
				pUnwrapper->DoUnwrap(true);
				RKProgress::Get().SetProgress(100);

				m_WorkingList[Index]->SetNewUVs(pUnwrapper);
				pUnwrapper->DeleteVertsAndFaces();
				delete pUnwrapper;
			}
		}
	}

	RKProgress::Get().CloseProgressWindow();
}



void RoadkillCmd::ClearWorkingList()
{
	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		WorkingMesh* pMesh = m_WorkingList[Index];
		delete pMesh;
	}
	m_WorkingList.clear();
}



/****************************************************
*		Set up the Synatax for Roadkill				*
****************************************************/

MSyntax RoadkillCmd::newSyntax()
{
	MSyntax syntax = MPxPolyTweakUVCommand::newSyntax();

	syntax.addFlag(kFunctionFlag, kFunctionFlagLong, MSyntax::kString); 
	syntax.addFlag(kScaleFlag, kScaleFlagLong, MSyntax::kUnsigned);
	syntax.addFlag(kUprightFlag, kUprightFlagLong, MSyntax::kUnsigned);
	syntax.addFlag(kInOneFlag, kInOneFlagLong, MSyntax::kUnsigned);
//	syntax.addFlag(kIslandMatchFlag, kIslandMatchFlagLong, MSyntax::kDouble);
//	syntax.addFlag(kIslandEdgeMatchFlag, kIslandEdgeMatchFlagLong, MSyntax::kDouble);

	return syntax;
}



/****************************************************
*		Parse the Synatax for Roadkill				*
****************************************************/

MStatus RoadkillCmd::parseSyntax (MArgDatabase &argData)
{
	MStatus status;

	MString Function;
	status = argData.getFlagArgument(kFunctionFlagLong, 0, Function);
	if(status != MS::kSuccess) return status;
	const char* funcType = Function.asChar();


	gMaintainUVs = false;
	gMatch = false;
	gEdgeLen = false;
	gSymMatch = false;

	int Result = 0;


	MString MaintainUVsExists = "optionVar -exists rkMaintainUVs ";
	MGlobal::executeCommand(MaintainUVsExists, Result);

	if (Result)
	{
		int MaintainUVs = false;
		MString getMaintainUVs = "optionVar -q rkMaintainUVs";
		MGlobal::executeCommand(getMaintainUVs, MaintainUVs);
		gMaintainUVs = MaintainUVs;
	}


	MString matchSymShellExists = "optionVar -exists rkIslandSymMatch ";
	MGlobal::executeCommand(matchSymShellExists, Result);

	if(Result)
	{
		int SymMatch = false;
		MString getMatchingSym = "optionVar -q rkIslandSymMatch";
		MGlobal::executeCommand(getMatchingSym, SymMatch);
		gSymMatch = SymMatch;
	}


	MString matchShellExists = "optionVar -exists rkIslandMatch";
	MGlobal::executeCommand(matchShellExists, Result);

	if(Result)
	{
		int Match = false;
		MString getMatching = "optionVar -q rkIslandMatch";
		MGlobal::executeCommand(getMatching, Match);
		gMatch = Match;
	}


	MString edgeMatchFlag = "optionVar -exists rkIslandEdgeMatch";
	MGlobal::executeCommand(edgeMatchFlag, Result);

	if(Result)
	{
		int EdgeMatch = false;
		MString getEdgeMatch = "optionVar -q rkIslandEdgeMatch";
		MGlobal::executeCommand(getEdgeMatch, EdgeMatch);
		gEdgeLen = EdgeMatch;
	}


	MString angleTolExists = "optionVar -exists rkIslandMatchAnglesTol";
	MGlobal::executeCommand(angleTolExists, Result);

	if(Result)
	{
		int AngleTol = 4;
		MString getAngleTol = "optionVar -q rkIslandMatchAnglesTol";
		MGlobal::executeCommand(getAngleTol, AngleTol);
		gMatchAngleTolerence = (double)AngleTol;
		if(gMatchAngleTolerence < 0.0) gMatchAngleTolerence = 0.0;
		gMatchAngleTolerence *= 0.0174532925; 
	}



	MString edgeTolExists = "optionVar -exists rkIslandMatchEdgeTol";
	MGlobal::executeCommand(edgeTolExists, Result);

	if(Result)
	{
		int EdgeTol = 4;
		MString getEdgeTol = "optionVar -q rkIslandMatchEdgeTol";
		MGlobal::executeCommand(getEdgeTol, EdgeTol);
		gMatchEdgeTolerence = (double)EdgeTol;
		if(gMatchEdgeTolerence < 0.0) gMatchEdgeTolerence = 0.0;
		gMatchEdgeTolerence = gMatchEdgeTolerence / 100.0;
	}



	if(strcmp(funcType, "GEOM") == 0) gFunction = 0;
	if(strcmp(funcType, "ORGA") == 0) gFunction = 1;
	if(strcmp(funcType, "SCALEPACK") == 0) gFunction = 2;
	if(strcmp(funcType, "PACK") == 0) gFunction = 3;
	if(strcmp(funcType, "MINIMISE") == 0) gFunction = 4;
	if(strcmp(funcType, "S2S") == 0) gFunction = 6;
	if(strcmp(funcType, "STRA") == 0) gFunction = 7;
	if(strcmp(funcType, "SYMCUT") == 0) gFunction = 8;


	int Scale = 10;
	int Upright = true;
	int AllInOne = false;

	status = argData.getFlagArgument(kScaleFlag, 0, Scale);
	if(status == MS::kSuccess) gScale = Scale;

	status = argData.getFlagArgument(kUprightFlag, 0, Upright);
	if(status == MS::kSuccess) gUpright = Upright;

	status = argData.getFlagArgument(kInOneFlag, 0, AllInOne);
	if(status == MS::kSuccess) gAllInOnePack = AllInOne;

	gOuterBorderSwitch = 0;
	gOuterBorderPixels = 4;
	gOuterBorderTMapSize = 1024;

	MGlobal::executeCommand("optionVar -exists rkOuterBorder", Result);
	if(Result) MGlobal::executeCommand("optionVar -q rkOuterBorder", gOuterBorderSwitch);

	MGlobal::executeCommand("optionVar -exists rkOuterBorderPixels", Result);
	if(Result) MGlobal::executeCommand("optionVar -q rkOuterBorderPixels", gOuterBorderPixels);

	MGlobal::executeCommand("optionVar -exists rkOuterBorderTMap", Result);
	if(Result) MGlobal::executeCommand("optionVar -q rkOuterBorderTMap", gOuterBorderTMapSize);

	Setup();

	return MS::kSuccess;
}




/************************************************
*	In comes a selected mesh, unwrap and set	*
*	its UVs										*
************************************************/

MStatus RoadkillCmd::getTweakedUVs(const MObject& meshObj, MIntArray& uvList, MFloatArray & uPos, MFloatArray & vPos)
{
	if(gSetupFailed == true) return MS::kFailure;
	
	MFnMesh mesh( meshObj );
	MFnDependencyNode depNodeFn;
	depNodeFn.setObject(meshObj);
	MString meshNodeName = depNodeFn.name();

	MDagPath dPath;
	mesh.getPath(dPath);
	MObject	meshTrans = dPath.transform();
	MFnTransform	trans(meshTrans);
	MString transName = trans.fullPathName();


	WorkingMesh* pThisMesh = findMesh(transName);
	if(pThisMesh == NULL) return MS::kFailure;

	if(gFunction != 8)
	{
		if(pThisMesh->m_SelectedEdges.length() != 0)
		{
			int NumWorkingUVs = pThisMesh->m_UCoords.length();
			uvList.setLength(NumWorkingUVs);
			uPos.setLength(NumWorkingUVs);
			vPos.setLength(NumWorkingUVs);

			for(int Index = 0; Index < NumWorkingUVs; Index++)
			{
				uvList[Index] = Index;
				uPos[Index] = pThisMesh->m_UCoords[Index];
				vPos[Index] = pThisMesh->m_VCoords[Index];
			}

		}
		else
		{
			// only work on selected UVs
			int NumWorkingUVs = uvList.length();
			uPos.setLength(NumWorkingUVs);
			vPos.setLength(NumWorkingUVs);

			for(int Index = 0; Index < NumWorkingUVs; Index++)
			{
				int UVIndex = uvList[Index];
				uPos[Index] = pThisMesh->m_UCoords[UVIndex];
				vPos[Index] = pThisMesh->m_VCoords[UVIndex];
			}
		}
	}

	removeMesh(meshNodeName);

	return MS::kSuccess;
}




WorkingMesh* RoadkillCmd::findMesh(MString& meshNodeName)
{
	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		const char* pName = meshNodeName.asChar();
		const char* pTestName = m_WorkingList[Index]->m_meshNodeName.asChar();

		if(m_WorkingList[Index]->m_meshNodeName == meshNodeName) return(m_WorkingList[Index]);
	}

	return(NULL);
}




void RoadkillCmd::removeMesh(MString& meshNodeName)
{
	vector<WorkingMesh*> newList;

	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		WorkingMesh* pMesh = m_WorkingList[Index];

		if(pMesh->m_meshNodeName == meshNodeName)
		{
			delete pMesh;
		}
		else
		{
			newList.push_back(pMesh);
		}
	}

	m_WorkingList.clear();
	m_WorkingList = newList;
}








void RoadkillCmd::SymmetryCut()
{
	MSelectionList StoreSelList;

	MGlobal::getActiveSelectionList(StoreSelList);

	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		m_WorkingList[Index]->SymmetryCut();				// find history!
	}

	MGlobal::setActiveSelectionList(StoreSelList);
}




/*******************************************************
*	Source is always the last mesh in the WorkingList  *
*******************************************************/

void RoadkillCmd::ScaleToSource()
{
	double SourceRatio  = 0.0;

	for(int Index = 0; Index < m_WorkingList.size(); Index++)
	{
		m_WorkingList[Index]->GetUVToFaceRatio();
		SourceRatio = m_WorkingList[Index]->m_UvToFaceRatio;
	}

	if(SourceRatio == 0.0)
	{
		displayError( "Roadkill: Unable to scale to source" );
	}

	for(int Index = 0; Index < m_WorkingList.size()-1; Index++)
	{
		double testScale = m_WorkingList[Index]->m_UvToFaceRatio;
		if(testScale != 0.0)
		{
			double NewScale = sqrt(SourceRatio / testScale);
			m_WorkingList[Index]->ScaleUVs(NewScale);
		}
	}
}
