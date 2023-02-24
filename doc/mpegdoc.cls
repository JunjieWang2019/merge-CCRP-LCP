%%
%% mpegdoc: Copyright (c) 2011 BBC Research & Development
%% All rights reserved.
%%
%% Permission to use, copy, modify, and/or distribute this software for any
%% purpose with or without fee is hereby granted, provided that the above
%% copyright notice and this permission notice appear in all copies.
%%
%% Neither the names of the BBC, ITU, ISO, IEC nor the names of its
%% contributors may be used to endorse or promote products derived from
%% this software without specific prior written permission.
%%
%% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
%% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
%% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
%% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
%% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
%% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
%% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mpegdoc}[2022/11/23 MPEG WG7 document class]

\LoadClassWithOptions{article}

% typeset using Times New Roman
\RequirePackage{ifxetex}
\ifxetex
  % use the correct fonts when using XeTeX
  \RequirePackage{fontspec}
  \defaultfontfeatures{Mapping=tex-text}
  \setmainfont{Times New Roman}
\else
  \RequirePackage{times}
\fi

% expand the margins
\RequirePackage{geometry}
\geometry{tmargin=1.6cm,lmargin=1in,rmargin=2.5cm,bmargin=1in,nohead}

% require graphicx for loading the logos in the header
\RequirePackage{graphicx}

% don't number the abstract
\renewenvironment{abstract}{\section*{Abstract}}{}

% format the date in iso style
\RequirePackage{datetime}
\renewcommand{\dateseparator}{-}
\newdateformat{@mpegdoc@yyyymmdddate}{%
	\THEYEAR\dateseparator\twodigit{\THEMONTH}\dateseparator\twodigit{\THEDAY}}

% all pages only have a footer
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\renewcommand{\headrulewidth}{0pt}
\fancyhead{}
\cfoot{\thepage}
\rfoot{Date saved: \@mpegdoc@yyyymmdddate\today}

% macros for document metadata
\RequirePackage[normalem]{ulem}
\def\@mpegdocnum{\uline{xxx}}
\def\@mpegdocstatus{\uline{TODO: Add \texttt{\textbackslash{}mpegdocstatus}}}
\def\@mpegdocpurpose{\uline{TODO: Add \texttt{\textbackslash{}mpegdocpurpose}}}
\def\@mpegdocsource{\uline{TODO: Add \texttt{\textbackslash{}mpegdocsource}}}
\def\@mpegdocserialnumber{\uline{TODO: Add \texttt{\textbackslash{}mpegdocserialnumber}}}
\def\@mpegmeeting{\uline{TODO: Add \texttt{\textbackslash{}mpegmeeting}}}
\newcommand{\mpegdocnum}[1]{\def\@mpegdocnum{\mbox{#1}}}
\newcommand{\mpegdocstatus}[1]{\def\@mpegdocstatus{\mbox{#1}}}
\newcommand{\mpegdocpurpose}[1]{\def\@mpegdocpurpose{\mbox{#1}}}
\newcommand{\mpegdocsource}[1]{\def\@mpegdocsource{\mbox{#1}}}
\newcommand{\mpegdocserialnumber}[1]{\def\@mpegdocserialnumber{\mbox{#1}}}
\newcommand{\mpegmeeting}[1]{\def\@mpegmeeting{\mbox{#1}}}

% affiliation block in the author list
\newenvironment{affiliation}{%
	\begin{tabular}[t]{@{}>{\renewcommand{\\}{, }}l@{}}%
}{%
	\end{tabular}%
}

\RequirePackage{array}
\RequirePackage{hyperref}

\RequirePackage{ifthen}
\newcounter{mpeg@author@column}
\newcommand{\@mpeg@switch@author@column}{%
	\ifthenelse{\value{mpeg@author@column} = 1}{%
		% do nothing
	}{%
		\addtocounter{mpeg@author@column}{1}%
		&
	}%
}

\newcommand*{\email}[1]{\@mpeg@switch@author@column\textit{#1}}
\newcommand*{\tel}[1]{\@mpeg@switch@author@column\textrm{#1}}

% count total number of pages including cover page)
\RequirePackage{totpages}

% Make a title at the top of the current page.
\renewcommand\maketitle{%
	\thispagestyle{empty}
	% none of the title is put in the header (it is too difficult to
	% make the page match the word style in this case)
	%\vspace*{-1.7cm}
	\newcommand{\@strutb}{\rule{0pt}{2.5ex}}

	\begin{raggedright}
	\rule{0pt}{1.5cm}% strut incase images aren't loaded
	\IfFileExists{logo-iso-iec.png}{\includegraphics[height=1.5cm]{logo-iso-iec.png}}{}\\
		\vspace*{-1.2cm}
		\hfill \textbf{\Large ISO/IEC JTC 1/SC 29/WG 7 \- }\textbf{\Huge\@mpegdocnum}\\[0ex]
		\vspace*{-.2cm}
		\hfill \rule{0.75\linewidth}{0.4pt}

		%\rule{\linewidth}{0.4pt}

	\vspace{\baselineskip}

	\begin{center}
		\Large
		\sffamily
		\begin{tabular*}{\linewidth}{|l@{\extracolsep{\fill}}c@{\extracolsep{\fill}}r|}
			\hline
			\@strutb & \textbf{ISO/IEC JTC 1/SC 29/WG 7} & \\
			\@strutb & \textbf{MPEG 3D Graphics Coding and Haptics} & \\
			\@strutb & \textbf{Convenorship: AFNOR (France)} & \\
			\hline
		\end{tabular*}
	\end{center}

	\vspace{0.5\baselineskip}

	\begin{tabular}{lp{0.78\textwidth}}
		\@strutb \textbf{Document Type:} & Output Document \\
		\\
		\@strutb \textbf{Title:} & \@title \\
		\\
		\@strutb \textbf{Status:} & \@mpegdocstatus \\
		\\
		\@strutb \textbf{Date of document:} & \@mpegdoc@yyyymmdddate\today \\
		\\
		\@strutb \textbf{Source:} & ISO/IEC JTC 1/SC 29/WG 7 \\
		\\
		\@strutb \textbf{Expected action:} & None \\
		\\
		\@strutb \textbf{Action due date:} & None \\
		\\
		\@strutb \textbf{No. of pages:} & \pageref{TotPages} \\
		\\
		\@strutb \textbf{Email of Convenor:} & \href{mailto:marius.preda@imt.fr}{marius.preda@imt.fr} \\
		\\
		\@strutb \textbf{Committee URL:} & {\small\url{https://isotc.iso.org/livelink/livelink/open/jtc1sc29wg7}} \\

	\end{tabular}

	%\rule{\linewidth}{0.4pt}

	\clearpage
	\setcounter{page}{1}

	\begin{center}
		\Large
		\textbf{INTERNATIONAL ORGANISATION FOR STANDARDIZATION}\\[0ex]
		\textbf{ORGANISATION INTERNATIONALE DE NORMALISATION}\\[0ex]
		\textbf{ISO/IEC JTC 1/SC 29/WG 7}\\[0ex]
		\textbf{MPEG 3D GRAPHICS CODING AND HAPTICS}\\[0ex]
		\-\\[0ex]
	\end{center}
		\hfill \textbf{\Large ISO/IEC JTC 1/SC 29/WG 7 \- }\textbf{\Huge\@mpegdocnum}\\[0ex]
		\hfill \textbf{\Large\@mpegmeeting}

	\vspace{0.5\baselineskip}

	\begin{tabular}{lp{0.78\textwidth}}
	\@strutb \textbf{Title} & \textbf{\@title} \\
	%\@strutb \it Purpose: & \@mpegdocpurpose \\
	%\@strutb \it Editors: & %
	%	\setcounter{mpeg@author@column}{0}
	%	\let\@and\\
	%	\renewcommand{\and}{\@and\setcounter{mpeg@author@column}{0}}
	%	\newcommand{\@NLtoPAR}{\renewcommand{\\}{\par}}
	%	\begin{tabular}[t]{@{}>{\@NLtoPAR}p{3in}|>{\@NLtoPAR}p{2in}@{}}%
	%		\@author
	%	\end{tabular}\\
	\@strutb \textbf{Source} & \textbf{\@mpegdocsource} \\
	\@strutb \textbf{Status} & \textbf{\@mpegdocstatus} \\
	\@strutb \textbf{Serial Number} & \textbf{\@mpegdocserialnumber} \\
	\end{tabular}
	\end{raggedright}

	% draw a short horizontal line to delimit the title from body
	{\center\rule{0.35\textwidth}{1pt}\\}
}

%
\RequirePackage{parskip}
\RequirePackage[compact]{titlesec}
%\titlespacing{\section}{0pt}{*0}{*0}
%\titlespacing{\subsection}{0pt}{*0}{*0}
%\titlespacing{\subsubsection}{0pt}{*0}{*0}

\titlespacing*{\section}
{0pt}{5ex}{2ex}
\titlespacing*{\subsection}
{0pt}{5ex}{2ex}
\titlespacing*{\subsection}
{0pt}{5ex}{2ex}
