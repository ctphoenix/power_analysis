\documentclass{beamer}

\usepackage{amsmath}
\usepackage{amssymb} % Often needed for math symbols
\usepackage{graphicx}

% --- Define Title Elements HERE (before \begin{document}) ---
\title{Proposal for a Multi-Arm Biomarker RCT}
\subtitle{Executive Summary}
\author{}
\date{\today}

% Choose a theme (optional, default is fine)
%\usetheme{default} 

\begin{document}

% --- Slide 1: Title ---
\begin{frame}
  \titlepage % Display the title page using info defined above
\end{frame}

% --- Slide 2: Motivation ---
\begin{frame}{Motivation}
  \begin{itemize}
    \item \textbf{Goal:} Design a multi-arm RCT to validate treatment-predictive biomarkers in psychiatry. \pause
    \item \textbf{(1) Efficiency:} Validating known biomarker candidates is often more statistically efficient than discovering them \\{\it de novo}. \pause
    \item \textbf{(2) Future Data:} Yields a rich dataset suitable for \emph{post-hoc} development of complex predictive models (which treatment for whom?).
  \end{itemize}
\end{frame}

% --- Slide 3: Model Setup ---
\begin{frame}{Model \& Assumptions: Setup}
  \begin{itemize}
    \item \textbf{Trial Structure:} \(K\) distinct experimental treatment arms. \pause
    \item \textbf{Key Simplification:} Powering based on \textbf{one treatment, one biomarker} assumption. \pause
      \begin{itemize}
          \item Each arm \(k\) is paired prospectively with a unique biomarker \(X_k\). \pause
          \item Aligns with literature reporting \& simplifies power calculations. \pause
          \item (Acknowledges reality: one biomarker might affect multiple treatments).
      \end{itemize}
  \end{itemize}
\end{frame}

% --- Slide 4: Linear Model ---
\begin{frame}{Model \& Assumptions: Linear Model}
  \begin{itemize}
    \item \textbf{Goal:} Power the trial to detect significance of parameters in: \pause
          \[ Y_i = \beta_0 + \sum_{k=1}^K [ \beta_{1k} T_{ik} + \beta_{2k} X_{ik} + \beta_{3k} (T_{ik} \cdot X_{ik}) ] + \text{error}_i \] \pause % Using display math here is fine in pure Beamer
    \item \(T_{ik}\) = 1 if patient \(i\) in arm \(k\), 0 otherwise. \pause
    \item \(X_{ik}\) = Biomarker value for patient \(i\) associated with arm \(k\).
  \end{itemize}
\end{frame}

% --- Slide 5: Coefficient Interpretation ---
\begin{frame}{Model \& Assumptions: Coefficient Interpretation}
  \begin{itemize}
    \item \(\beta_0\) : Average outcome for the reference treatment group. \pause
    \item \(\beta_{1k}\) : Main effect of treatment \(k\) vs. reference (traditional RCT target). \pause
    \item \(\beta_{2k}\) : General association of biomarker \(X_k\) with outcome (prognostic effect; nuisance/confounder control). \pause
    \item \(\mathbf{\beta_{3k}}\) : \textbf{Interaction effect} (Primary Target!) - How biomarker \(X_k\) modulates treatment \(k\)'s effect. Goal is accurate, significant estimation.
  \end{itemize}
\end{frame}

% --- Slide 6: Power Basis ---
\begin{frame}{Power Estimates: Basis}
  \begin{itemize}
    \item \textbf{Source:} Simulations based on plausible effect sizes (\(\beta_{3k}\)) derived from literature. \pause
    \item \textbf{Examples Modeled (K=4):} \pause
      \begin{itemize}
          \item iAPF (EEG) for rTMS (e.g., \(\beta_3 \approx 0.2-0.3\)) \pause
          \item AUD History (EHR) for Ketamine (e.g., \(\beta_3 \approx 0.6\) or higher from OR) \pause
          \item Inflammatory Markers (Blood) for ECT (e.g., \(\beta_3 \approx 0.2\)) \pause
          \item Speech Latency (Voice) for Novel Agents (e.g., \(\beta_3 \approx 0.3\)) \pause
      \end{itemize}
    \item \textbf{Caveats:} Literature estimates are rough (differing populations, disease criteria, outcome measures, reported stats). Use with caution.
  \end{itemize}
\end{frame}

% --- Slide 7: Power Definitions ---
\begin{frame}{Power Estimates: Definitions & Control}
  \begin{itemize}
    \item \textbf{Statistical Power:} Probability of correctly detecting a true effect (interaction \(\beta_{3k} \neq 0\)) when it exists. Target often 80\%. \pause
    \item \textbf{Family-Wise Error Rate (FWER):} Probability of making \emph{at least one} false positive discovery (claiming \(\beta_{3k} \neq 0\) when it's truly 0) across all K tests. \pause 
    \item \textbf{Control Method:} Used Holm's procedure (conservative) to control FWER at 5\%.
  \end{itemize}
\end{frame}

% --- Slide 8: Simulation Results (Figure Only) ---
\begin{frame}{Power Estimates: Simulation Results}
  \begin{figure}
      \centering
      \includegraphics[width=\textwidth]{figs/biomarker_power_holm.png}
      \caption{Power per biomarker vs. sample size \emph{per arm} (K=4, Holm FWER control).}
  \end{figure} 
\end{frame}

% --- Slide 9: Simulation Results (Text Discussion) ---
\begin{frame}{Power Estimates: Simulation Results}
  \begin{itemize}
    \item \textbf{Observation 1:} Required sample size varies greatly by interaction strength (\(\beta_{3k}\)). \pause
    \item \textbf{Example (80\% Power):} \pause
      \begin{itemize}
          \item Strong (\(\beta_{3k} \approx 0.6\)): ~50 subjects/arm \pause
          \item Weak (\(\beta_{3k} \approx 0.2\)): ~350 subjects/arm \pause
      \end{itemize}
    \item \textbf{Observation 2:} If goal is to validate \emph{all} K biomarkers, trial must be sized for the \emph{weakest} interaction (350/arm -> Total N = {\bf 1400}).
  \end{itemize}
\end{frame}

% --- Slide 10: Considerations ---
\begin{frame}{Considerations & Future Directions}
  \begin{itemize}
    \item \textbf{Biomarker Correlation:} Assumed biomarkers are correlated (via latent factor U). \pause
      \begin{itemize}
          \item Higher correlation generally \emph{decreases} power, potentially requiring larger sample sizes than shown if correlation is strong. (Estimating this pre-trial is valuable). \pause
      \end{itemize}
    \item \textbf{Adaptive Designs?:} Could potentially \emph{reduce} sample size (e.g., drop futile arms/biomarkers early, enrich promising ones). \pause
      \begin{itemize}
          \item Trade-off: Increased operational and statistical complexity. \pause
      \end{itemize}
    \item \textbf{Forthcoming: Predictive Modeling Simulations:} Planning extended simulations focusing on building \emph{predictive models} using the full feature set (many biomarkers, instruments) to generate patient-specific treatment profiles (see separate plan).
  \end{itemize}
\end{frame}

% --- Slide 11: End ---
\begin{frame}
  \centering
  \Huge End
\end{frame}

\end{document}