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

    // strategy tab
    float pre_test_infect_prob{1.};
    int mode;
    std::string mode_str;
    int time_passed;
    int quarantine;
    std::vector<int> t_test{};
    QString test_type;

    // parameters tab
    std::vector<float> residence_times_mean;
    std::vector<float> residence_times_lev;
    std::vector<float> residence_times_uev;

    float fraction_asymtomatic;
    float sensitivity;
    float specificity;

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

    void collect_data(MainWindow *parent);
    std::vector<int> collect_t_test(std::vector<QCheckBox*> boxes);

    Eigen::VectorXf calc_rates(std::vector<float> times,
                               std::vector<int> comp);
    Eigen::MatrixXf calc_S(int n);
    Eigen::MatrixXf calc_A(Eigen::MatrixXf S_,
                           Eigen::VectorXf r);
    Eigen::MatrixXf calc_X(int delay,
                           int qrntn,
                           Eigen::MatrixXf A_,
                           Eigen::VectorXf states);
    Eigen::MatrixXf assemble_phases(Eigen::MatrixXf X_,
                                    std::vector<int> comp);

    QtCharts::QChartView* create_plot();
    QTableWidget* create_table();

    void create_result_log();
    void write_row_result_log(QTableWidget*);
    void update_result_log();

    float calculate_strategy_result(Eigen::MatrixXf matrix);
signals:

public slots:
};

#endif // SIMULATION_H
