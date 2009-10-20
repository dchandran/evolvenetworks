/****************************************************************************

Copyright (c) 2008 Deepak Chandran
Contact: Deepak Chandran (dchandran1@gmail.com)
See COPYRIGHT.TXT

This file defines the class that is used to create new threads in the 
Tinkercell main window. The threads can be associated with a dialog that provides
users with the option to terminate the thread.


****************************************************************************/
#include "CThread.h"
#include <QVBoxLayout>
#include <QDockWidget>
#include <QDir>
#include <QLabel>
#include <QSemaphore>
#include <QCoreApplication>
#include <QtDebug>

namespace Tinkercell
{
	/******************
	PROCESS THREAD
	*******************/

	QWidget * ProcessThread::dialog(QWidget * mainWindow, ProcessThread * newThread, const QString& text, QIcon icon)
	{
		QWidget * dialog = new QDialog(mainWindow);

		dialog->hide();

		dialog->move(mainWindow->pos() + QPoint(10,10));
		dialog->setWindowIcon(icon);

		QHBoxLayout * layout = new QHBoxLayout;
		QPushButton * killButton = new QPushButton("Terminate Program");
		connect(killButton,SIGNAL(released()),dialog,SLOT(accept()));
		QLabel * label = new QLabel(text + tr(" is Running..."));

		layout->addWidget(label);
		layout->addWidget(killButton);
		dialog->setWindowTitle(tr("Program Running"));

		dialog->setLayout(layout);
		dialog->setWindowFlags(Qt::Dialog);
		dialog->setAttribute(Qt::WA_DeleteOnClose,true);

		connect(killButton,SIGNAL(released()),newThread,SLOT(terminate()));
		connect(newThread,SIGNAL(finished()),dialog,SLOT(close()));
		connect(newThread,SIGNAL(started()),dialog,SLOT(show()));

		return dialog;
	}

	ProcessThread::ProcessThread(const QString& exe, const QString& args, QWidget * main)
		: QThread(main)
	{
		this->args = args;
		this->exe = exe;
		
		connect(this,SIGNAL(terminated()),this,SLOT(stopProcess()));
		connect(this,SIGNAL(finished()),this,SLOT(stopProcess()));
	}

	void ProcessThread::run()
	{
		if (!exe.isEmpty())
		{
			//setPriority(QThread::LowestPriority);
			connect(this,SIGNAL(terminated()),&process,SLOT(kill()));
			process.start(exe,QStringList() << args);
			process.waitForFinished();
			errStream = process.readAllStandardError();
			outputStream = process.readAllStandardOutput();
			//ConsoleWindow::error(errors);
		}
	}

	void ProcessThread::stopProcess()
	{
		if (process.state() != QProcess::NotRunning)
			process.close();
	}

	ProcessThread::~ProcessThread()
	{
		stopProcess();
	}
	
	QString ProcessThread::output() const
	{
		return outputStream;
	}
	
	QString ProcessThread::errors() const
	{
		return errStream;
	}

}
