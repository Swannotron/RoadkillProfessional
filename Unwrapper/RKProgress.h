#ifndef RKPROGRESS
#define RKPROGRESS


class RKProgress
{
public:
	~RKProgress() {};

	static RKProgress& Get();

	void NewProgressWindow();
	void CloseProgressWindow();
	void SetText(const char* pNewText);
	void SetProgress(int Amount);

private:
	RKProgress() {};
};




#endif

