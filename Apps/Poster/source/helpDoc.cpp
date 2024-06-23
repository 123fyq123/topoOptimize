#include "helpDoc.h"
#include "ui_helpDoc.h"
#include <iostream>


HelpDoc::HelpDoc(QWidget* parent) :
	QDialog(parent),
	ui(new Ui::HelpDoc)
{

	ui->setupUi(this);

}

HelpDoc::~HelpDoc()
{
	std::cout << "[DEBUG] " << "进行了析构\n";
	delete ui;
}
