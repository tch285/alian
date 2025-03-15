#include "apprun.hh"
#include <TGApplication.h> // Added include for TGApplication

namespace alian
{

	AppRun::AppRun(const char *name, int *argc, char **argv)
			: TApplication(name, argc, argv),
				fMain(nullptr),
				fQuitButton(nullptr)
	{
		// Create the main window with specified dimensions.
		fMain = new TGMainFrame(gClient->GetRoot(), 500, 200);

		// Create a button labeled "Quit". When clicked, it will call the Quit slot.
		fQuitButton = new TGTextButton(fMain, "Quit");
		// Set layout to stretch the button horizontally.
		fQuitButton->SetTextJustify(36); // Center text horizontally and vertically.
		fQuitButton->Connect("Clicked()", "alian::AppRun", this, "Quit()");
		// Add the button to the window with centering layout hints.
		fMain->AddFrame(fQuitButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 50, 50, 50, 50));
		fMain->SetWindowName("Quit Application");

		// Map subwindows, resize, and display the main window.
		fMain->MapSubwindows();
		fMain->Resize();
		fMain->MapWindow();
	}

	AppRun::~AppRun()
	{
		fMain->Cleanup();
		delete fMain;
	}

	int AppRun::run()
	{
		TApplication::Run(); // Run() returns void, so we simply call it.
		return 0;						 // Return 0 after exiting the event loop.
	}

	void AppRun::Quit()
	{
		fMain->CloseWindow();
		TApplication::Terminate(0);
	}

	void wait_for_quit()
	{
		int argc = 0;
		char **argv = nullptr;
		// Create a TGApplication to initialize ROOT's GUI.
		// Note: TGApplication derives from TApplication.
		TGApplication app("TGApp", &argc, argv);
		// Create our application window.
		AppRun myApp("AppRun", &argc, argv);
		myApp.run();
	}

} // namespace alian