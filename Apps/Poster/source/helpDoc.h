#pragma once

#include <QtWidgets/QMainWindow>
#include <QDialog>


namespace Ui {
	class HelpDoc;
}

class HelpDoc : public QDialog
{
	Q_OBJECT

public:
	explicit HelpDoc(QWidget* parent = nullptr);
	~HelpDoc();

private:
	Ui::HelpDoc* ui;
};