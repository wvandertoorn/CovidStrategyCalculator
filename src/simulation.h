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
    explicit Simulation();
    explicit Simulation(MainWindow *parent);
    void run();

protected:
    MainWindow* m_parent = nullptr;
    static std::vector<int> sub_compartments;
    static int nr_compartments;

    bool time_passed_known{false};
    int time_passed{0};
    int quarantine{11};
    std::vector<int> t_test{};

    // parameters
    std::vector<float> residence_times_mean = {2., 3.5, 4.81, 5};
    std::vector<float> residence_times_lev = {1., 1.5, 2.81, 2};
    std::vector<float> residence_times_uev = {5., 7.5, 10.81, 8};
    float pre_test_infect_prob{1.};
    float percentage_asymt{.2};
    float pcr_sensitivity{.70};
    float pcr_specificity{.995};

    // simulation
    Eigen::VectorXf initial_states;
    Eigen::MatrixXf result_matrix_mean;
    Eigen::MatrixXf result_matrix_lev;
    Eigen::MatrixXf result_matrix_uev;

    void collect_data(MainWindow *parent);
    std::vector<int> collect_t_test(std::vector<QCheckBox*> boxes);
    Eigen::VectorXf calc_rates(std::vector<float> times,
                               std::vector<int> comp);
    Eigen::MatrixXf calc_S(int n);
    Eigen::MatrixXf calc_A(Eigen::MatrixXf S_,
                           Eigen::VectorXf r);
    Eigen::MatrixXf calc_X(float delay,
                           float qrntn,
                           Eigen::MatrixXf A_,
                           Eigen::VectorXf states);
    Eigen::MatrixXf assemble_phases(Eigen::MatrixXf X_,
                                    std::vector<int> comp);

    QtCharts::QChartView* create_plot();
    QTableWidget* create_table();
    float calculate_strategy_result(Eigen::MatrixXf matrix);
    void create_result_log();
    void write_row_result_log(QTableWidget*);
    void update_result_log();
    void output_results();

signals:

public slots:
};

#endif // SIMULATION_H
