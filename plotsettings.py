from	astropy.io import ascii
import	matplotlib
from	matplotlib.path import Path
import	matplotlib.patches as patches
import	matplotlib.path as mpath
import	matplotlib.patheffects as PathEffects
import	matplotlib.pylab as plt
from	matplotlib import rc, rcParams
import	numpy as np
import	platform

left, width                     = 0, 1
bottom, height                  = 0, 1

right                           = left + width
top                             = bottom + height

fontsize                        = 32
label_size                      = fontsize - 2
legend_size                     = fontsize - 12

rcParams['xtick.major.width'] 	= 1.25
rcParams['ytick.major.width'] 	= 1.25

rcParams['xtick.minor.width'] 	= 1.25
rcParams['ytick.minor.width']	= 1.25

rcParams['xtick.major.size'] 	= 8
rcParams['ytick.major.size'] 	= 8

rcParams['xtick.minor.size'] 	= 5
rcParams['ytick.minor.size'] 	= 5

rcParams['xtick.labelsize'] 	= label_size
rcParams['ytick.labelsize'] 	= label_size

rcParams['legend.numpoints'] 	= 1
rcParams['legend.scatterpoints']= 1

rcParams['savefig.dpi']         = '1200'
rcParams['savefig.bbox']        = 'tight'

rcParams['xtick.direction']		= 'in'
rcParams['ytick.direction']		= 'in'

rcParams['xtick.top']		    = True
rcParams['xtick.bottom']		= True
rcParams['ytick.left']		    = True
rcParams['ytick.right'] 		= True

#rcParams['font.family'] 		= 'serif'
#rcParams['mathtext.fontset'] 	= 'dejavuserif'
rcParams['text.usetex'] 		= True
rc('text.latex', preamble='\\usepackage{cmbright}')


rcParams['axes.linewidth']		= 2
rcParams['axes.labelsize']		= fontsize

box 				            = dict(facecolor='white', pad=5, linewidth=0.0, alpha=0.75)

if len(platform.mac_ver()[0]) > 0:
	matplotlib.use('MacOSX')
else:
	matplotlib.use('Qt5Agg')

# Colours
# Solarize colour scheme

color_back                      = '#fdf6e3'
color_green_light				= '#c2cc7f'
color_green						= '#859900'
color_green_dark				= '#424c00'
color_base03					= '#002b36'
color_base02					= '#073642'
color_base01					= '#586e75'
color_base00					= '#657b83'
color_base0						= '#839496'
color_base1						= '#93a1a1'
color_base2						= '#eee8d5'
color_base3						= '#fdf6e3'
color_yellow_light				= '#dac47f'
color_yellow					= '#b58900'
color_orange					= '#cb4b16'
color_red  						= '#dc322f'
color_magenta					= '#d33682'
color_violet					= '#6c71c4'
color_blue  					= '#268bd2'
color_cyan  					= '#2aa198'
color_green 					= '#859900'
	
# Vigit colour scheme

vigit_color_1                   = '#0D3D56'
vigit_color_2                   = '#0F5B78'
vigit_color_3                   = '#117899'
vigit_color_4                   = '#1496BB'
vigit_color_5                   = '#5CA794'
vigit_color_6                   = '#A3B86C'
vigit_color_7                   = '#EBC944'
vigit_color_8                   = '#EDAA38'
vigit_color_9                   = '#F08C2D'
vigit_color_10                  = '#F26D21'
vigit_color_11                  = '#D94E20'
vigit_color_12                  = '#C02F1D'
vigit_grey_1                    = '#131516'
vigit_grey_2                    = '#373D3F'
vigit_grey_3                    = '#555F61'
vigit_grey_4                    = '#707C80'
vigit_grey_5                    = '#8C979A'
vigit_grey_6                    = '#A7B0B2'
vigit_grey_7                    = '#C1C7C9'
vigit_grey_8                    = '#DADEDF'
vigit_grey_9                    = '#F2F3F4'

vigit_color_1_light             = '#869eaa'
vigit_color_2_light             = '#87adbb'
vigit_color_3_light             = '#88bbcc'
vigit_color_4_light             = '#89cadd'
vigit_color_5_light             = '#add3c9'
vigit_color_6_light             = '#d1dbb5'
vigit_color_7_light             = '#f5e4a1'
vigit_color_8_light             = '#f6d49b'
vigit_color_9_light             = '#f7c596'
vigit_color_10_light            = '#f8b690'
vigit_color_11_light            = '#eca68f'
vigit_color_12_light            = '#df978e'
vigit_color_6_dark              = '#515c36'
vigit_color_9_dark              = '#784616'
vigit_color_10_dark             = '#793610'
vigit_color_11_dark             = '#6c2710'
vigit_color_12_dark             = '#60170e'


colors_vigit					= [vigit_color_1, vigit_color_2, vigit_color_3, vigit_color_4, vigit_color_5,
								vigit_color_6, vigit_color_7, vigit_color_8, vigit_color_9, vigit_color_10,
								vigit_color_11, vigit_color_12, vigit_grey_1, vigit_grey_2, vigit_grey_3,
								vigit_grey_4, vigit_grey_5, vigit_grey_6, vigit_grey_7, vigit_grey_8, vigit_grey_9
								]