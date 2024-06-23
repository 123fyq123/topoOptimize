#pragma once

#include <QtWidgets/QMainWindow>
#include <QDialog>


namespace Ui {
    class AboutDialog;
}

class AboutDialog : public QDialog
{
	Q_OBJECT

public:
	explicit AboutDialog(QWidget* parent = nullptr);
	~AboutDialog();

private:
	Ui::AboutDialog* ui;
};