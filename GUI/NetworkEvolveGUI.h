#ifndef NETWORK_EVOLVE_GUI_H
#define NETWORK_EVOLVE_GUI_H

#include <QApplication>
#include <QCoreApplication>
#include <QWidget>
#include <QMainWindow>
#include <QSettings>
#include <QLabel>
#include <QTextEdit>
#include <QPushButton>
#include <QSplitter>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QGroupBox>
#include <QComboBox>
#include <QLineEdit>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QToolBar>
#include <QFile>
#include <QProcess>
#include <QTextStream>
#include <QFileInfo>
#include <QDir>
#include <QFileDialog>
#include "CodeEditor.h"
#include "SyntaxHighlighter.h"

namespace NetworkEvolutionLib
{
	class MainWindow : public QMainWindow
	{
		Q_OBJECT
		
	public:
		MainWindow();
		~MainWindow();
		QSize sizeHint() const;
		
	private:
	
		QComboBox * fitnessComboBox;
		Tinkercell::CodeEditor * codeEditor;
		QWidget * setupNetworkOptions();
		void setupMassActionNetwork(QTreeWidget*);
		void setupEnzymeNetwork(QTreeWidget*);
		void setupProteinNetwork(QTreeWidget*);
		void setupGeneticNetwork(QTreeWidget*);
		
		QWidget * setupGAOptions();
		QWidget * setupEditor();
		QString init();
		QString mainFunction();
		QString callbackFunction();
		QComboBox* comboBox();
		
		/*params to fill in*/
		
		double uni_uni, uni_bi, bi_uni, bi_bi, no_reactant, no_product;
		double ma_min_constant, ma_max_constant;
		double ma_mutate_constants, ma_mutate_remove_reaction, ma_mutate_add_reaction;
		
		double prot_min_ka, prot_min_vmax, prot_min_total;
		double prot_max_ka, prot_max_vmax, prot_max_total;
		double prot_mutate_rewire, prot_mutate_parameter, prot_mutate_total, prot_mutate_addremove;
		
		double enzyme_min_kcat, enzyme_min_log_keq, enzyme_min_alpha, enzyme_min_h, enzyme_min_s_half, enzyme_min_p_half;
		double enzyme_max_kcat, enzyme_max_log_keq, enzyme_max_alpha, enzyme_max_h, enzyme_max_s_half, enzyme_max_p_half;
		double 	enzyme_mutate_enzyme, 
				enzyme_mutate_k_cat, 
				enzyme_mutate_k_eq, 
				enzyme_mutate_alpha,
				enzyme_mutate_h, 
				enzyme_mutate_s_half, 
				enzyme_mutate_p_half, 
				enzyme_mutate_remove, 
				enzyme_mutate_add;
	
		double grn_init_inflow, grn_init_cost_per_protein;
		int grn_min_complex_size;
		int grn_max_complex_size;
		double grn_min_Ka, grn_min_Vmax, grn_min_degradation;
		double grn_max_Ka, grn_max_Vmax, grn_max_degradation;
		double grn_mutate_Ka, grn_mutate_Vmax, grn_mutate_complex, grn_mutate_add_gene, grn_mutate_remove_gene;

		int species_min, species_max, reactions_min, reactions_max;
		double crossover_rate;
		double init_iv;
		double mutate_iv;
		bool lineageTracking;
		
		double mass_action_prob, enzyme_prob, protein_net_prob, grn_prob;
		
		QString codeFile, logFile;
		int runs, generations, popSz;
		QString seeds;
		QString compileCommand;
		
		int bestNetworkFitness1, bestNetworkScript1,
			bestNetworkSize1, bestNetworkLineage1,
			allFitness1, allNetworkLineage1;
			
		int bestNetworkFitness2, bestNetworkScript2,
			bestNetworkSize2, bestNetworkLineage2,
			allFitness2, allNetworkLineage2, seeds2;

		double max_fitness;
		
	private slots:
	
		void clear();
		void save();
		void run();
		void reset();
		void fitnessSelected(const QString&);
		
		void useMassAction(bool use) { mass_action_prob = 1.0 * (int)use; }
		void useEnzyme(bool use) { enzyme_prob = 1.0 * (int)use; }
		void useProtein(bool use) { protein_net_prob = 1.0 * (int)use; }
		void useGRN(bool use) { grn_prob = 1.0 * (int)use; }
		
		void setMassActionProb(double value) { mass_action_prob = value; }
		void setEnzymeProb(double value) { enzyme_prob = value; }
		void setProteinProb(double value) { protein_net_prob = value; }
		void setGRNProb(double value) { grn_prob = value; }		
		
		void setUniUni(double value) { uni_uni = value; }
		void setUniBi(double value) { uni_bi = value; }
		void setBiUni(double value) { bi_uni = value; }
		void setBiBi(double value) { bi_bi = value; }
		void setNoReactant(double value) { no_reactant = value; }
		void setNoProduct(double value) { no_product = value; }
		
		void set_ma_min_constant(double value) { ma_min_constant = value; }
		void set_ma_max_constant(double value) { ma_max_constant = value; }
		void set_ma_mutate_constants(double value) { ma_mutate_constants = value; }
		void set_ma_mutate_remove_reaction(double value) { ma_mutate_remove_reaction = value; }
		void set_ma_mutate_add_reaction(double value) { ma_mutate_add_reaction = value; }
		
		void set_enzyme_min_kcat(double value) { enzyme_min_kcat = value; }
		void set_enzyme_min_log_keq(double value) { enzyme_min_log_keq = value; }
		void set_enzyme_min_alpha(double value) { enzyme_min_alpha = value; }
		void set_enzyme_min_h(double value) { enzyme_min_h = value; }
		void set_enzyme_min_s_half(double value) { enzyme_min_s_half = value; }
		void set_enzyme_min_p_half(double value) { enzyme_min_p_half = value; }
		
		void set_enzyme_max_kcat(double value) { enzyme_max_kcat = value; }
		void set_enzyme_max_log_keq(double value) { enzyme_max_log_keq = value; }
		void set_enzyme_max_alpha(double value) { enzyme_max_alpha = value; }
		void set_enzyme_max_h(double value) { enzyme_max_h = value; }
		void set_enzyme_max_s_half(double value) { enzyme_max_s_half = value; }
		void set_enzyme_max_p_half(double value) { enzyme_max_p_half = value; }
		
		void set_enzyme_mutate_enzyme(double value) { enzyme_mutate_enzyme = value; }
		void set_enzyme_mutate_k_cat(double value) { enzyme_mutate_k_cat = value; }
		void set_enzyme_mutate_k_eq(double value) { enzyme_mutate_k_eq = value; }
		void set_enzyme_mutate_alpha(double value) { enzyme_mutate_alpha = value; }
		void set_enzyme_mutate_h(double value) { enzyme_mutate_h = value; }
		void set_enzyme_mutate_s_half(double value) { enzyme_mutate_s_half = value; }
		void set_enzyme_mutate_p_half(double value) { enzyme_mutate_p_half = value; }
		void set_enzyme_mutate_remove(double value) { enzyme_mutate_remove = value; }
		void set_enzyme_mutate_add(double value) { enzyme_mutate_add = value; }
		
		void set_prot_min_ka(double value) { prot_min_ka = value; }
		void set_prot_min_vmax(double value) { prot_min_vmax = value; }
		void set_prot_min_total(double value) { prot_min_total = value; }
		void set_prot_max_ka(double value) { prot_max_ka = value; }
		void set_prot_max_vmax(double value) { prot_max_vmax = value; }
		void set_prot_max_total(double value) { prot_max_total = value; }
		
		void set_prot_mutate_rewire(double value) { prot_mutate_rewire = value; }
		void set_prot_mutate_parameter(double value) { prot_mutate_parameter = value; }
		void set_prot_mutate_total(double value) { prot_mutate_total = value; }
		void set_prot_mutate_addremove(double value) { prot_mutate_addremove = value; }
		
		void set_grn_init_inflow(double value) { grn_init_inflow = value; }
		void set_grn_init_cost_per_protein(double value) { grn_init_cost_per_protein = value; }
		void set_grn_min_complex_size(int value) { grn_min_complex_size = value; }
		void set_grn_max_complex_size(int value) { grn_max_complex_size = value; }
		
		void set_grn_min_Ka(double value) { grn_min_Ka = value; }
		void set_grn_min_Vmax(double value) { grn_min_Vmax = value; }
		void set_grn_min_degradation(double value) { grn_min_degradation = value; }
		void set_grn_max_Ka(double value) { grn_max_Ka = value; }
		void set_grn_max_Vmax(double value) { grn_max_Vmax = value; }
		void set_grn_max_degradation(double value) { grn_max_degradation = value; }
		
		void set_grn_mutate_Ka(double value) { grn_mutate_Ka = value; }
		void set_grn_mutate_Vmax(double value) { grn_mutate_Vmax = value; }
		void set_grn_mutate_complex(double value) { grn_mutate_complex = value; }
		void set_grn_mutate_add_gene(double value) { grn_mutate_add_gene = value; }
		void set_grn_mutate_remove_gene(double value) { grn_mutate_remove_gene = value; }
		
		void setCodeFile(const QString& file) 
		{ 
			compileCommand.replace(codeFile,file);
			codeFile = file;
		}
		void setLogFile(const QString& file) { logFile = file; }
		void setSeed(const QString& s) { seeds = s; }
		void setRuns(int i) { runs = i; }
		void setGenerations(int i) { generations = i; }
		void setPopSz(int i) { popSz = i; }
		
		void set_min_species(int value) { species_min = value; }
		void set_min_reactions(int value) { reactions_min = value; }
		void set_max_species(int value) { species_max = value; }
		void set_max_reactions(int value) { reactions_max = value; }
		void set_crossover_rate(double value) { crossover_rate = value; }
		void set_init_iv(double value) { init_iv = value; }
		void set_mutate_iv(double value) { mutate_iv = value; }
		void set_lineageTracking(bool value) { lineageTracking = value; }
		
		void setCompileCommand(const QString& s) { compileCommand = s; }
		
		void setBestNetworkFitness1(bool value) { bestNetworkFitness1 = (int)value; }
		void setBestNetworkScript1(bool value) { bestNetworkScript1 = (int)value; }
		void setBestNetworkSize1(bool value) { bestNetworkSize1 = (int)value; }
		void setBestNetworkLineage1(bool value) { bestNetworkLineage1 = (int)value; }
		void setAllFitness1(bool value) { allFitness1 = (int)value; }
		void setAllNetworkLineage1(bool value) { allNetworkLineage1 = (int)value; }
		
		void setBestNetworkFitness2(bool value) { bestNetworkFitness2 = (int)value; }
		void setBestNetworkScript2(bool value) { bestNetworkScript2 = (int)value; }
		void setBestNetworkSize2(bool value) { bestNetworkSize2 = (int)value; }
		void setBestNetworkLineage2(bool value) { bestNetworkLineage2 = (int)value; }
		void setAllFitness2(bool value) { allFitness2 = (int)value; }
		void setAllNetworkLineage2(bool value) { allNetworkLineage2 = (int)value; }
		
		void showSeed(bool value) { seeds = (int)value; }
		void setMaxFitness(double value) { max_fitness = (int)value; }
	};
}

#endif

