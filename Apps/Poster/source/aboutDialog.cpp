#include "aboutDialog.h"
#include "ui_aboutDialog.h"
#include <iostream>


AboutDialog::AboutDialog(QWidget* parent) :
	QDialog(parent),
	ui(new Ui::AboutDialog)
{

	ui->setupUi(this);
	
}

AboutDialog::~AboutDialog()
{
	std::cout << "[DEBUG] " << "进行了析构\n";
	delete ui;
}
