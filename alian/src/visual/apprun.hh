#ifndef ALIAN_APPRUN_HH
#define ALIAN_APPRUN_HH

#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>

namespace alian
{

	class AppRun : public TApplication
	{
	public:
		AppRun(const char *name, int *argc, char **argv);
		virtual ~AppRun();
		int run();
		void Quit(); // Slot for button click

	private:
		TGMainFrame *fMain;
		TGTextButton *fQuitButton;
	};

	void wait_for_quit();

} // namespace alian

#endif // ALIAN_APPRUN_HH