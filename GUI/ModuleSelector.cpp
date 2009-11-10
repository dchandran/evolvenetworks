#include <QApplication>
#include <QCoreApplication>
#include "ModuleSelector.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
	
	QString appDir = QCoreApplication::applicationDirPath();
	
    /*QFile styleFile(appDir + QString("/networkevolve.qss"));
	
	if (styleFile.open(QFile::ReadOnly | QFile::Text))
    {
        app.setStyleSheet(styleFile.readAll());
        styleFile.close();
    }*/

    NetworkEvolutionLib::ModuleSelector widget;

    widget.show();

	return app.exec();
}

namespace NetworkEvolutionLib
{
	ModuleSelector::ModuleSelector(QWidget * parent): QWidget(parent)
	{
		QTabWidget * tabWidget = new QTabWidget;
		
		QString appDir = QCoreApplication::applicationDirPath();
		
		QDir dir(appDir);
		if (dir.cd("modules"))
		{
			QFileInfoList fileInfoList = dir.entryInfoList();
			
			for (int i=0; i < fileInfoList.size(); ++i)
			{
				if (fileInfoList[i].isDir() && !fileInfoList[i].baseName().isEmpty())
				{
					if (dir.cd(fileInfoList[i].baseName()))
					{
						QFileInfoList fileInfoList2 = dir.entryInfoList(QStringList() << "*.png" << "*.PNG");
						
						tabWidget->addTab(
							generateWidget(fileInfoList2),
							fileInfoList[i].baseName().replace("_"," ")
						);
						
						dir.cdUp();
					}
				}
			}			
		}
		
		QHBoxLayout * layout = new QHBoxLayout;
		layout->addWidget(tabWidget);
		layout->setContentsMargins(0,0,0,0);
		setLayout(layout);
	}
	
	QWidget* ModuleSelector::generateWidget(const QFileInfoList& files)
	{
		QWidget * widget = new QWidget;
		QGridLayout * layout = new QGridLayout;
		
		for (int i=0, j=0, k=0; k < files.size(); ++j, ++k)
		{
			QToolButton * button = new QToolButton;
			QCheckBox * checkbox = new QCheckBox;
			QVBoxLayout * vlayout = new QVBoxLayout;
			vlayout->addWidget(button);
			vlayout->addWidget(checkbox);
			button->setIcon(QIcon(files[k].absoluteFilePath()));
			//button->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
			button->setIconSize(QSize(200,200));
			checkbox->setText(files[k].baseName().replace(tr("_"),tr(" ")));
			button->setStyleSheet(tr("QToolButton:pressed { background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,stop: 0 #dadbde, stop: 1 #f6f7fa); }"));
			if (j > 3)
			{
				++i;
				j = 0;
			}
			connect(button,SIGNAL(clicked()),checkbox,SLOT(toggle()));
			layout->addLayout(vlayout,i,j);
		}
		
		widget->setLayout(layout);
		return widget;
	}
}
