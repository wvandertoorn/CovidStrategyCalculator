#ifndef SIMULATION_H
#define SIMULATION_H

#include "mainwindow.h"

#include <Eigen/Dense>

#include <QCheckBox>
#include <QMainWindow>
#include <QTableWidget>

#include <QtCharts/QChartView>

class MainWindow;

class Simulation : public QObject
{
    Q_OBJECT
public:
    explicit Simulation(MainWindow *parent);
    explicit Simulation(MainWindow *parent, Eigen::MatrixXf init_states, float prerisk);
    void run();
    std::tuple<std::vector<Eigen::MatrixXf>,
               std::vector<Eigen::MatrixXf>> run_prevalence(std::vector<float> week_incidences);
    void output_results();

protected:
    // model structure
    MainWindow* m_parent = nullptr;
    static std::vector<int> sub_compartments;
    static int nr_compartments;
    static float t_inf;


    // strategy tab
    float pre_test_infect_prob{1.};
    int mode;
    std::string mode_str;
    int time_passed;
    int quarantine;
    bool use_symptomatic_screening;
    std::vector<int> t_test{};
    QString test_type;

    // parameters tab
    std::vector<float> residence_times_mean;
    std::vector<float> residence_times_lev;
    std::vector<float> residence_times_uev;

    float sensitivity;
    float specificity;

    float fraction_asymtomatic;

    // simulation
    Eigen::VectorXf initial_states;
    Eigen::MatrixXf X_mean;
    Eigen::MatrixXf X_lev;
    Eigen::MatrixXf X_uev;
    Eigen::MatrixXf result_matrix_mean;
    Eigen::MatrixXf result_matrix_lev;
    Eigen::MatrixXf result_matrix_uev;
    Eigen::MatrixXf assay_detectibility_worst_case;
    Eigen::MatrixXf assay_detectibility_best_case;
    Eigen::MatrixXf assay_detectibility_mean_case;

    float risk_T_postprocedure_mean, risk_T_postprocedure_lev, risk_T_postprocedure_uev;
    float risk_T_preprocedure_mean, risk_T_preprocedure_lev, risk_T_preprocedure_uev;
    float fold_RR_mean, fold_RR_lev, fold_RR_uev;
    std::vector<int> time_for_plot;

    void collect_data(MainWindow *parent);
    std::vector<int> collect_t_test(std::vector<QCheckBox*> boxes);

    Eigen::VectorXf calc_rates(std::vector<float> times,
                               std::vector<int> comp);
    Eigen::VectorXf FOR_vector(std::vector<int> comp);
    Eigen::MatrixXf calc_S(int n);
    Eigen::MatrixXf calc_A(Eigen::MatrixXf S_,
                           Eigen::VectorXf r,
                           std::vector<int> comp);
    Eigen::MatrixXf calc_X(int delay,
                           int qrntn,
                           Eigen::MatrixXf A_,
                           Eigen::VectorXf states);
    float calc_risk_at_T(Eigen::MatrixXf A,
                        Eigen::VectorXf states,
                        float time_T);
    Eigen::MatrixXf assemble_phases(Eigen::MatrixXf X_,
                                    std::vector<int> comp);

    Eigen::MatrixXf risk_node_to_relative_residual_risk(Eigen::MatrixXf risk, float risk_T);
    QtCharts::QChartView* create_plot(Eigen::MatrixXf mean,
                                      Eigen::MatrixXf uev,
                                      Eigen::MatrixXf lev,
                                      std::vector<int> time_range_for_plot);
    QTableWidget* create_table();

    void create_result_log();
    void write_row_result_log(QTableWidget*);
    void update_result_log();

    float calculate_strategy_result(Eigen::MatrixXf matrix);
    Eigen::MatrixXf calculate_assay_sensitivity();

    std::tuple<Eigen::MatrixXf, std::vector<int>, float> calculate_strategy_with_test(Eigen::MatrixXf A);
signals:

public slots:
};

#endif // SIMULATION_H
