#include <maya/MPxPolyTweakUVCommand.h>
#include <maya/MFnComponentListData.h>

#include <vector>
#include <WorkingMesh.h>

using namespace std;



//class MSyntax;
//class MArgDatabase;
//class MIntArray;
//class MFloatArray;
//class MObject;


#define kFunctionFlag				"-f" 
#define kFunctionFlagLong			"-function" 
#define kScaleFlag					"-s" 
#define kScaleFlagLong				"-scale" 
#define kUprightFlag				"-u" 
#define kUprightFlagLong			"-upright" 
#define kInOneFlag					"-i"
#define kInOneFlagLong				"-inone"
//#define kIslandMatchFlag			"-m" 
//#define kIslandMatchFlagLong		"-match"
//#define kIslandEdgeMatchFlag		"-e" 
//#define kIslandEdgeMatchFlagLong	"-edge" 






class RoadkillCmd : public MPxPolyTweakUVCommand
{
public:
	RoadkillCmd();
	~RoadkillCmd();

	static void *creator();
	static MSyntax newSyntax();

	bool		isUndoable() const;


	MStatus parseSyntax(MArgDatabase &argData);

//	MStatus getTweakedUVs(const MFnMesh & mesh,	MIntArray & uvList,	MFloatArray & uPos, MFloatArray & vPos);

	MStatus getTweakedUVs(const MObject & mesh,	MIntArray & uvList,	MFloatArray & uPos,	MFloatArray & vPos);

	void Setup();

private:

//	void UnwrapMesh(WorkingMesh* pWorkingMesh);
	void ScaleToSource();
	void SymmetryCut();

	void ClearWorkingList();
	WorkingMesh* findMesh(MString& meshNodeName);
	void removeMesh(MString& meshNodeName);

	vector<WorkingMesh*> m_WorkingList;
};
