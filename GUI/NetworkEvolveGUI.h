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
	
		QWidget * setupNetworkOptions();
		void setupMassActionNetwork(QTreeWidget*);
		void setupEnzymeNetwork(QTreeWidget*);
		void setupProteinNetwork(QTreeWidget*);
		void setupGeneticNetwork(QTreeWidget*);
		
		QWidget * setupGAOptions();
		QWidget * setupEditor();
		QString init();
		
		/*params to fill in*/
		
		double uni_uni, uni_bi, bi_uni, bi_bi, no_reactant, no_product;
		double ma_init_max_constant;
		double ma_mutate_constants, ma_mutate_remove_reaction, ma_mutate_add_reaction;
		
		double prot_init_ka, prot_init_vmax, prot_init_total;
		double prot_mutate_rewire, prot_mutate_parameter, prot_mutate_total, prot_mutate_addremove;
		
		double enzyme_init_max_kcat, enzyme_init_max_log_keq, enzyme_init_max_alpha, enzyme_init_max_h, enzyme_init_max_s_half, enzyme_init_max_p_half;
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
		int grn_init_max_complex_size, grn_init_Ka, grn_init_Vmax, grn_init_degradation;
		double grn_mutate_Ka, grn_mutate_Vmax, grn_mutate_complex, grn_mutate_add_gene, grn_mutate_remove_gene;

		int species, reactions;
		double crossover_rate;
		double init_iv;
		double mutate_iv;
		bool lineageTracking;
		
		double mass_action_prob, enzyme_prob, protein_net_prob, grn_prob;
		
		/*
		void setDistributionOfMassActionNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product);
		void setRateConstantForMassActionNetwork(double max_rate_constant);
		void setMutationRatesForMassActionNetwork(double prob_mutate_constants, double prob_mutate_remove_reaction, double prob_mutate_add_reaction);
		
		void setRateConstantsForProteinInteractionNetwork(double ka, double vmax, double total);
		void setMutationRatesForProteinInteractionNetwork(double rewire, double parameter, double total, double addremove);
		
		void setRateConstantsForEnzymeNetwork(double max_kcat, double max_log_keq, double max_alpha, double max_h, double max_s_half, double max_p_half);
		void setMutationRatesForEnzymeNetwork(double enzyme, 
									  double k_cat, 
									  double k_eq, 
									  double alpha,
									  double h, 
									  double s_half, 
									  double p_half, 
									  double remove, 
									  double add);
	
		void setResourceRestriction(double inflow, double cost_per_protein);
		void setRateConstantsForGeneRegulationNetwork(int max_complex_size, double Ka, double Vmax, double degradation);
		void setMutationRatesForGeneRegulationNetwork(double Ka, double Vmax, double complex, double add_gene, double remove_gene);

		setInitialNetworkSize(int,int);
		void setCrossoverRate(double);
		void setAverageInitialValue(double);
		void setMutationRateOfInitialValues(double);
		void lineageTrackingON();
		void lineageTrackingOFF();
		void setNetworkTypeProbability(int, double);
		void setNetworkType(int);
		
		*/
	private slots:
	
		void run();
		void reset();
		void quit();
		
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
		
		void set_ma_init_max_constant(double value) { ma_init_max_constant = value; }
		void set_ma_mutate_constants(double value) { ma_mutate_constants = value; }
		void set_ma_mutate_remove_reaction(double value) { ma_mutate_remove_reaction = value; }
		void set_ma_mutate_add_reaction(double value) { ma_mutate_add_reaction = value; }
		
		void set_enzyme_init_max_kcat(double value) { enzyme_init_max_kcat = value; }
		void set_enzyme_init_max_log_keq(double value) { enzyme_init_max_log_keq = value; }
		void set_enzyme_init_max_alpha(double value) { enzyme_init_max_alpha = value; }
		void set_enzyme_init_max_h(double value) { enzyme_init_max_h = value; }
		void set_enzyme_init_max_s_half(double value) { enzyme_init_max_s_half = value; }
		void set_enzyme_init_max_p_half(double value) { enzyme_init_max_p_half = value; }
		
		void set_enzyme_mutate_enzyme(double value) { enzyme_mutate_enzyme = value; }
		void set_enzyme_mutate_k_cat(double value) { enzyme_mutate_k_cat = value; }
		void set_enzyme_mutate_k_eq(double value) { enzyme_mutate_k_eq = value; }
		void set_enzyme_mutate_alpha(double value) { enzyme_mutate_alpha = value; }
		void set_enzyme_mutate_h(double value) { enzyme_mutate_h = value; }
		void set_enzyme_mutate_s_half(double value) { enzyme_mutate_s_half = value; }
		void set_enzyme_mutate_p_half(double value) { enzyme_mutate_p_half = value; }
		void set_enzyme_mutate_remove(double value) { enzyme_mutate_remove = value; }
		void set_enzyme_mutate_add(double value) { enzyme_mutate_add = value; }
		
		
	};
}

#endif

