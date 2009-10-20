/****************************************************************************

Copyright (c) 2008 Deepak Chandran
Contact: Deepak Chandran (dchandran1@gmail.com)
See COPYRIGHT.TXT

This file defines the class that is used to create new threads in the 
Tinkercell main window. The threads can be associated with a dialog that provides
users with the option to terminate the thread.


****************************************************************************/

#ifndef TINKERCELL_CTHREAD_H
#define TINKERCELL_CTHREAD_H

#include <QMainWindow>
#include <QTextEdit>
#include <QSyntaxHighlighter>
#include <QHash>
#include <QTextCharFormat>
#include <QDialog>
#include <QCompleter>
#include <QListWidget>
#include <QThread>
#include <QToolBar>
#include <QTabWidget>
#include <QTableWidget>
#include <QComboBox>
#include <QPushButton>
#include <QActionGroup>
#include <QLibrary>
#include <QProcess>
#include <QProgressBar>
#include <QItemDelegate>

#ifdef Q_WS_WIN
#define MY_EXPORT __declspec(dllexport)
#else
#define MY_EXPORT
#endif

namespace Tinkercell
{
	/*! \brief This class is used to run a process (command + args) as a separate thread as a separate thread
	\ingroup core
	*/
	class MY_EXPORT ProcessThread : public QThread
	{
		Q_OBJECT
	public:
		/*! \brief constructor -- used to initialize the main window, the command name and the args for the command
		* \param QString command
		* \param QString arguments
		* \param QWidget main window
		*/
		ProcessThread(const QString&, const QString& ,QWidget * main);
		/*! \brief get the results (output stream) from the process
		* \return QString output*/
		virtual QString output() const;
		/*! \brief get the errors (error stream) from the process
		* \return QString output*/
		virtual QString errors() const;
		/*! \brief destructor -- free the library that this thread loaded*/
		virtual ~ProcessThread();
		/*! \brief  creates a dialog that shows the name of the running thread and a button for terminating the thread
		* \param QWidget main window
		* \param ProcessThread
		* \param QString text to display
		* \param QIcon icon to display
		*/
		static QWidget * dialog(QWidget *, ProcessThread*, const QString& text = QString("Process"), QIcon icon = QIcon());
	protected slots:
		/*! \brief unload the library (if loaded) and delete it*/
		virtual void stopProcess();
	protected:
		/*! \brief the name of the executable*/
		QString exe;
		/*! \brief the arguments*/
		QString args;
		/*! \brief the output from the process*/
		QString outputStream;
		/*! \brief the error from the process*/
		QString errStream;
		/*! \brief Tinkercell's main window*/
		QProcess process;
		/*! \brief initializes the function pointers through the main window and then runs the target function*/
		virtual void run();

	};

}

#endif
