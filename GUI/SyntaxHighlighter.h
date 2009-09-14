#ifndef C_SYNTAX_HIGHLIGHTER_H
#define C_SYNTAX_HIGHLIGHTER_H

#include <QMainWindow>
#include <QTextEdit>
#include <QHash>
#include <QTextCharFormat>
#include <QDialog>
#include <QCompleter>
#include <QListWidget>
#include <QThread>
#include <QToolBar>
#include <QActionGroup>
#include <QRegExp>
#include <QSyntaxHighlighter>

namespace NetworkEvolutionLib
{
	class CSyntaxHighlighter : public QSyntaxHighlighter
    {
		 Q_OBJECT

	 public:
		 CSyntaxHighlighter(QTextDocument *parent = 0);

	 protected:
		 void highlightBlock(const QString &text);

	 private:
		 struct HighlightingRule
		 {
			 QRegExp pattern;
			 QTextCharFormat format;
		 };
		 QVector<HighlightingRule> highlightingRules;

		 QRegExp commentStartExpression;
		 QRegExp commentEndExpression;

		 QTextCharFormat keywordFormat;
		 QTextCharFormat classFormat;
		 QTextCharFormat singleLineCommentFormat;
		 QTextCharFormat multiLineCommentFormat;
		 QTextCharFormat quotationFormat;
		 QTextCharFormat functionFormat;
		 QTextCharFormat loopFormat1, loopFormat2;
	 };
}
#endif
