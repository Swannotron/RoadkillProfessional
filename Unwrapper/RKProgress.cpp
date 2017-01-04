#include <maya/MProgressWindow.h>
#include <maya/MString.h>


#include "RKProgress.h"



RKProgress& RKProgress::Get()
{
	static RKProgress	singleton;
	return singleton;
}


void RKProgress::NewProgressWindow()
{
	MProgressWindow::reserve();
	MProgressWindow::setTitle("Roadkill");
	MProgressWindow::setProgressRange(0, 100);
	MProgressWindow::setProgress(0);
	MProgressWindow::startProgress();
}


void RKProgress::CloseProgressWindow()
{
	MProgressWindow::endProgress();
}


void RKProgress::SetText(char* pNewText)
{
	MProgressWindow::setProgressStatus(pNewText);
}


void RKProgress::SetProgress(int Amount)
{
	MProgressWindow::setProgress(Amount);
}
