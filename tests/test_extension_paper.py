from graphpathology.extension_paper import clean_prose


def test_prose_spacing_preserves_latex_file_paths_ids_and_citations():
 raw=r'All48 fits: \includegraphics[width=\textwidth]{figures/ExtensionMatrix_fixed_20pct_image.pdf} \label{fig:fixed_20pct} \ref{fig:fixed_20pct} \cite{rampasek2022} TENX95 float64 ResNet18'
 result=clean_prose(raw)
 assert result.startswith('All 48 fits:')
 for token in [r'\includegraphics[width=\textwidth]{figures/ExtensionMatrix_fixed_20pct_image.pdf}',r'\label{fig:fixed_20pct}',r'\ref{fig:fixed_20pct}',r'\cite{rampasek2022}','TENX95','float64','ResNet18']:assert token in result
