import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class plot_diagrams:

    @classmethod
    def plot(self, hom, title):
        colors = ['C0', 'C1', 'C2']
        labels=['$H_0$','$H_1$', '$H_2$']

        fig = plt.figure(figsize=(6,6))
        gs  = fig.add_gridspec(nrows=2, ncols=2,
                            width_ratios=(3,1), height_ratios=(1,3),
                            wspace=0, hspace=0)
        ax_main = fig.add_subplot(gs[1,0])
        ax_main.set_ylabel("Death")
        ax_main.set_xlabel("Birth")

        ax_hisx = fig.add_subplot(gs[0,0],sharex=ax_main)
        ax_hisx.axis(False)

        ax_hisy = fig.add_subplot(gs[1,1],sharey=ax_main)
        ax_hisy.axis(False)

        for i in range(0, 3):
            x = hom[i][:, 0]
            y = hom[i][:, 1]
            ax_main.scatter(x, y, s=10, label=labels[i])

            if i != 0:
                ax_hisx.hist(x, bins=5, alpha=0.5, density=True, color=colors[i])

            ax_hisy.hist(y, bins=5, alpha=0.5, density=True, orientation='horizontal', color=colors[i])

        fig.suptitle(title)
        fig.legend()
        plt.show()


    @classmethod
    def plot_barcode_diagrams(self, betti_numbers, hom, title):
        fig = plt.figure(figsize=(6,12))
        gs = fig.add_gridspec(nrows=3,hspace=0)
        axs = gs.subplots(sharex=True)

        colors=['C0', 'C1', 'C2']
        labels=['$H_0$','$H_1$', '$H_2$']

        for i in range(0, betti_numbers.shape[0]):
            axs[i].set_ylabel(labels[i])
            axs[i].get_yaxis().set_ticks([])
            axs[i].label_outer()

            len_h = round(betti_numbers[i])
            for j in range(0, len_h):
                axs[i].plot(np.linspace(hom[i][j][0], hom[i][j][1], 2), [j, j], c=colors[i], linewidth=0.5)

        axs[0].set_title(title)
        axs[-1].set_xlabel("genetic distance")
        fig.legend()
        plt.show()

    
    @classmethod
    def add_trend_line_with_ci(self, ax, x, y, color, alpha=0.3):
        # Calculate the mean and standard deviation
        mean_y = np.mean(y, axis=0)
        std_y = np.std(y, axis=0)
        ci = 1.96 * (std_y / np.sqrt(len(y)))  # 95% confidence interval

        # Fit a polynomial of degree 1 (linear trend) to the data
        coefficients = np.polyfit(x, mean_y, 1)
        polynomial = np.poly1d(coefficients)
        trend_line = polynomial(x)

        # Plot the confidence interval
        ax.fill_between(x, trend_line - ci, trend_line + ci, color=color, alpha=alpha)


    @classmethod
    def plot_lines(self, ax, params, topoquant, xlab, ylab, main, lt):
        colors = {
            0: '#013366',  # H0
            1: '#ea7334',  # H1
            2: '#008b8b'   # H2
        }

        for hom in range(params["MAXDIM"] + 1):  # H0, H1, H2
            all_sims = []

            for sim_idx in range(params["SIMULATIONS"]):
                y = [topoquant[i][sim_idx][hom] for i in range(len(params[lt]))]

                ax.plot(params[lt], y, color=colors[hom],
                        alpha=0.4, linestyle='-', marker='o')
                all_sims.append(y)
            
            # Calculate the mean and confidence interval for the trend line
            self.add_trend_line_with_ci(ax, params[lt], np.array(all_sims), colors[hom])
            
            # Add a legend entry for each homology dimension
            ax.plot([], [], color=colors[hom], label=f'$H_{hom}$')
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(main)


    @classmethod
    def plot_topoquant_diagram(self, params, bn, blm, blv, lt, save_as_pdf, pdf):
        fig, axs = plt.subplots(1, 3, figsize=(21, 7)) 
        main_lab = 'Betti Numbers vs. '
        mbcl_lab = 'Barcode Length Mean vs. '
        vbcl_lab = 'Barcode Length Variance vs. '

        if lt == "SPALIST":
            plot_xlab = 'Sampling Ratio'
            main_lab += "Sampling Ratio"
            mbcl_lab += "Sampling Ratio"
            vbcl_lab += "Sampling Ratio"
        else:
            plot_xlab = 'Error Variance'
            main_lab += "Error Variance"
            mbcl_lab += "Error Variance"
            vbcl_lab += "Error Variance"

        self.plot_lines(axs[0], params, bn, xlab=plot_xlab,
                        ylab='Betti Numbers', main=main_lab, lt=lt)
        self.plot_lines(axs[1], params, blm, xlab=plot_xlab,
                        ylab='Mean Barcode Lengths', main=mbcl_lab, lt=lt)
        self.plot_lines(axs[2], params, blv, xlab=plot_xlab,
                        ylab='Barcode Lengths Variance', main=vbcl_lab, lt=lt)

        # Add a main title for the entire figure
        if lt == "SPALIST":
            fig.suptitle('Topological Quantities at Varying Sampling Ratios', fontsize=16)
        else:
            fig.suptitle('Topological Quantities at Varying Errors', fontsize=16)

        # Create a single legend for the entire figure
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', ncol=3)
        
        # Adjust layout to make room for the title
        plt.tight_layout(rect=[0, 0, 1, 0.96])  
        if save_as_pdf:
            pdf.savefig()
            plt.close()

        else:
            # For Jupyter Notebook handling
            plt.show()