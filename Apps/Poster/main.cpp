#include "iGamePoster.h"
#include <QtWidgets/QApplication>
#include <QtGUI/Qsurfaceformat.h>
#include <QtGUI/QPixmap>
#include <QtWidgets/Qsplashscreen.h>
#include <QtCore/Qtextcodec.h>

int main(int argc, char *argv[])
{
	QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
	//QCoreApplication::setAttribute (Qt::AA_UseDesktopOpenGL);
	QApplication a(argc, argv);
	a.processEvents ();
	QSurfaceFormat format;
	format.setRenderableType (QSurfaceFormat::OpenGL);
	/*format.setSamples (16);*/
	format.setDepthBufferSize (32);
	format.setStencilBufferSize (16);
	format.setVersion (4,1);
	//format.setProfile (QSurfaceFormat::CompatibilityProfile);
	format.setProfile (QSurfaceFormat::CoreProfile);
	QSurfaceFormat::setDefaultFormat (format);
	QPixmap pixmap ("./Resources/iGameLogo.png");
	QSplashScreen splash (pixmap);
	splash.show ();
	QTextCodec *codec = QTextCodec::codecForName ("GBK");
	iGamePoster w;
	w.setWindowTitle (codec->toUnicode ("iGame"));
    splash.finish (&w);
	w.show ();
	return a.exec();
}
