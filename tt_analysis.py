# Timeout-Tree Analysis
# Calculates optimal security delay (inactive lifetime) for a given timeout-tree (TT) with given input parameters, and calculates efficiency metrics given that security delay

# Usage:
# python3 tt_analysis.py -n YY
# Reads input parameters from file in_tt_analysisYY.csv
#     Row 1 (fixed parameter names): Ac,Ro,AS,MS
#     Row 2 (fixed parameter values): Ac_value,Ro_value,AS_value,MS_value
#     Row 3 (dynamic parameter names): Fe,Ex,Pr,Le,Va,Co
#     Row 4 (dynamic parameter values): Fe_1,Ex_1,Pr_1,Le_1,Va_1,Co_1
#     Row 5 (dynamic parameter values): Fe_2,Ex_2,Pr_2,Le_2,Va_2,Co_2
#     Row 6 (dynamic parameter values): Fe_3,Ex_3,Pr_3,Le_3,Va_3,Co_3
# Writes results to output file out_tt_analysisYY.csv
#     Row 1 (fixed parameter names): Ac,Ro,AS,MS
#     Row 2 (fixed parameter values): Ac_value,Ro_value,AS_value,MS_value
#     Row 3 (dynamic parameter and output metric names): Fe,Ex,Pr,Le,Va,Co,FractionTTLeaves,SecurityDelayBlocks,SecurityDelayYears,CapitalCost,CapitalEfficiency,OnchainFee,OnchainFeeFraction,ExpectedOnchainFee,ExpectedOverheadFraction
#     Row 4 (dynamic parameter values): Fe_1,Ex_1,Pr_1,Le_1,Va_1,Co_1,FractionTTLeaves_1,SecurityDelayBlocks_1,Security_DelayYears_1,CapitalCost_1,CapitalEfficiency_1,OnchainFee_1,OnchainFeeFraction_1,ExpectedOnchainFee_1,ExpectedOverheadFraction_1
#     Row 5 (dynamic parameter values): Fe_2,Ex_2,Pr_2,Le_2,Va_2,Co_2,FractionTTLeaves_2,SecurityDelayBlocks_2,Security_DelayYears_2,CapitalCost_2,CapitalEfficiency_2,OnchainFee_2,OnchainFeeFraction_2,ExpectedOnchainFee_2,ExpectedOverheadFraction_2
#     Row 6 (dynamic parameter values): Fe_3,Ex_3,Pr_3,Le_3,Va_3,Co_3,FractionTTLeaves_3,SecurityDelayBlocks_3,Security_DelayYears_3,CapitalCost_3,CapitalEfficiency_3,OnchainFee_3,OnchainFeeFraction_3,ExpectedOnchainFee_3,ExpectedOverheadFraction_3
# Provides analysis assuming feerates are of the form Fe/(1-x)^Ex where::
#     x is fraction of block used by TT leaves
#     Fe is feerate without any TT leaves
#     Ex is exponent controlling rate of growth of feerates as TT leaves are added to blocks

# Fixed Parameters: (used for all analyses) (specified in row 2 of input file)
# Ac: length (in blocks) of active lifetime of each TT
# Ro: length (in blocks) of rollover period of each TT (provides for casual user's unavailability)
# AS: average size (in vbytes) of transactions required to put one TT leaf onchain when all leaves in TT are put onchain
# MS: maximum size (in vbytes) of transactions required to put one TT leaf onchain when only one leaf in TT is put onchain

# Dynamic Parameters: (vary per analysis) (specified in rows 4+ of input file)
# Fe: feerate (in sats/vbyte) assuming 0% of block contains TT leaves
# Ex: exponent controlling rate of growth of feerates as TT leaves are added to blocks
# Pr: probability given TT is put on-chain
# Le: number of leaves of all TTs put together
# Va: value (in bitcoins) of all TT leaves put together, where each TT leaf has equal funds and during its active lifetime each leaf's funds are equally divided between:
#       1) casual user's immediate bitcoin,
#       2) casual user's Lightning balance, and
#       3) dedicated user's Lightning balance
# Co: cost of capital (in fraction of funder's capital/year) for allocating TT funder's capital to TT
#    Note: If TT leaves use hierarchical channels, funder is able to use their capital to route payments unrelated to the leaf's casual user, thus reducing the cost of capital

# e(x) = Expected cost per leaf for putting all Le leaves of all TTs onchain as a function of x is:
#    Va*SATOSHIS_PER_BITCOIN*Co*AS/(BLOCKS_PER_YEAR*BLOCKSIZE*x) + Pr*AS*calc_feerate(Fe, Ex, x)
# Expected cost is minimized by using binary search to minimize e(x) based on its derivative e'(x) where e'(x) is:
#    -Va*SATOSHIS_PER_BITCOIN*Co*AS/(BLOCKS_PER_YEAR*BLOCKSIZE*(x**2)) + Pr*AS*calc_feerate_derivative(Fe, Ex, x)

import sys
import math
import argparse
import csv

# Constants:
BLOCKSIZE = 4000000					# blocksize (in vbytes)
BLOCKS_PER_YEAR = 144*365.25				# number of blocks per year (on average)
SATOSHIS_PER_BITCOIN = 100000000			# number of satoshis per bitcoin
LOG_PRECISION = 50					# parameter that gives the number of binary search iterations in calculating the optimal value of x, where x is the fraction of each block used for putting TT leaves onchain

def calc_feerate(Fe, Ex, x):				# x is fraction of blocks used by TT leaves
    return(Fe / (1.0 - x)**Ex)

def calc_feerate_derivative(Fe, Ex, x):			# derivative of feerate with respect to x, where x is fraction of blocks used by TT leaves
    return(Fe * Ex / (1.0 - x)**(Ex+1))

def analyze_tt(static_params, dynamic_params, out_csv):
    Ac = int(static_params[0])
    Ro = int(static_params[1])
    AS = int(static_params[2])
    MS = int(static_params[3])
    Fe = float(dynamic_params[0])
    Ex = float(dynamic_params[1])
    Pr = float(dynamic_params[2])
    Le = int(dynamic_params[3])
    Va = int(dynamic_params[4])
    Co = float(dynamic_params[5])
    assert Ac >= 0.0
    assert Ro >= 0.0
    assert AS >= 0.0
    assert MS >= AS
    assert Fe >= 0.0
    assert Ex > 0.0
    assert 0.0 <= Pr <= 1.0
    assert Le >= 0.0
    assert Va >= 0.0
    assert 0.0 <= Co <= 1.0
    cuf = 2.0 * Va * SATOSHIS_PER_BITCOIN / (3.0 * Le)	# casual user's funds per leaf
    assert cuf > MS*Fe					# casual user's funds per leaf must be more than the maximum fee when fees are not increased due to congestion from TT leaves
    low_x = 0.0						# lower bound on allowable x value
    high_x = 1.0					# upper bound on allowable x value
    for i in range(LOG_PRECISION):			# first, determine largest x value for which max onchain fees do not exceed casual user's funds per leaf
        x = (low_x + high_x) / 2.0
        max_fee = MS*calc_feerate(Fe, Ex, x)
        if (max_fee > cuf):
            high_x = x
        else:
            low_x = x
    low_x = 0.0
    high_x = x
    for i in range(LOG_PRECISION):			# next, determine value of x that minimizes expected cost, subject to constraint that max onchain fees do not exceed casual user's funds per leaf
        x = (low_x + high_x) / 2.0
        exp_cost = Va*SATOSHIS_PER_BITCOIN*Co*AS/(BLOCKS_PER_YEAR*BLOCKSIZE*x) + Pr*AS*calc_feerate(Fe, Ex, x)
        exp_cost_deriv = -Va*SATOSHIS_PER_BITCOIN*Co*AS/(BLOCKS_PER_YEAR*BLOCKSIZE*(x**2)) + Pr*AS*calc_feerate_derivative(Fe, Ex, x)
        if (exp_cost_deriv > 0.0):
            high_x = x
        else:
            low_x = x
    max_fee = MS*calc_feerate(Fe, Ex, x)
    onchain_fee_fraction = max_fee / cuf
    ave_fee = AS*calc_feerate(Fe, Ex, x)
    expected_fee = Pr*ave_fee
    security_delay_blocks = int(math.ceil(Le*AS/(BLOCKSIZE*x)))
    security_delay_years = float(security_delay_blocks) / BLOCKS_PER_YEAR
    capital_cost = Va*SATOSHIS_PER_BITCOIN*Co*security_delay_blocks/(BLOCKS_PER_YEAR*Le)
    expected_cost = capital_cost + expected_fee
    capital_efficiency = (2.0/3.0) * Ac / (Ac + Ro + security_delay_blocks)
    expected_overhead_fraction = (capital_cost + expected_fee) / cuf
    out_csv.writerow([Fe,Ex,Pr,Le,Va,Co,x,security_delay_blocks,security_delay_years,capital_cost,capital_efficiency,max_fee,onchain_fee_fraction,expected_fee,expected_overhead_fraction])

# Main program
parser = argparse.ArgumentParser(description='Analyze timeout-tree security delays and efficiency metrics')
parser.add_argument('-n', dest='number', action='store', default='00', help='number portion of i/o file names')
args = parser.parse_args()
in_file_name = 'in_tt_analysis' + args.number + '.csv'
out_file_name = 'out_tt_analysis' + args.number + '.csv'
with open(in_file_name) as in_f:
    in_csv = csv.reader(in_f)
    static_headers = next(in_csv)
    expected_static_headers = ['Ac','Ro','AS','MS']
    assert expected_static_headers == static_headers
    static_params = next(in_csv)
    dynamic_headers = next(in_csv)
    expected_dynamic_headers = ['Fe','Ex','Pr','Le','Va','Co']
    assert expected_dynamic_headers == dynamic_headers
    with open(out_file_name, 'w') as out_f:
        out_csv = csv.writer(out_f)
        out_csv.writerow(['Timeout-Tree (TT) analysis'])
        out_csv.writerow([''])
        out_csv.writerow(static_headers)
        out_csv.writerow(['active','rollover','ave size','max size'])
        out_csv.writerow(['period','period','of txs','of txs'])
        out_csv.writerow(['(blocks)','(blocks)','(vbytes)','(vbytes)'])
        out_csv.writerow(static_params)
        out_csv.writerow([''])
        out_csv.writerow(['Fe','Ex','Pr','Le','Va','Co','FractionTTLeaves','SecurityDelayBlocks','SecurityDelayYears','CapitalCost','CapitalEfficiency','OnchainFee','OnchainFeeFraction','ExpectedOnchainFee','ExpectedOverheadFraction'])
        out_csv.writerow(['feerate','feerate','prob','leaves','value of all','cost of','fraction of block','delay for putting','delay for putting','capital cost','fraction of funder\'s','max fee','max fee per','expected fee per','capital cost plus'])
        out_csv.writerow(['base','exponent','TT put','across all TTs','leaves put together','capital','space devoted to','leaves onchain','leaves onchain','per leaf','funds used by','per onchain leaf','onchain leaf as fraction','leaf','expected fee per leaf as'])
        out_csv.writerow(['(sats/vbyte)','','onchain','','(BTC)','','leaves','(blocks)','(years)','(sats)','casual user','(sats)','of casual user\'s funds','(sats)','fraction of casual user\'s funds'])
        for dynamic_params in in_csv:
            analyze_tt(static_params, dynamic_params, out_csv)
        out_csv.writerow(['Fe','Ex','Pr','Le','Va','Co','FractionTTLeaves','SecurityDelayBlocks','SecurityDelayYears','CapitalCost','CapitalEfficiency','OnchainFee','OnchainFeeFraction','ExpectedOnchainFee','ExpectedOverheadFraction'])
        out_csv.writerow(['feerate','feerate','prob','leaves','value of all','cost of','fraction of block','delay for putting','delay for putting','capital cost','fraction of funder\'s','max fee','max fee per','expected fee per','capital cost plus'])
        out_csv.writerow(['base','exponent','TT put','across all TTs','leaves put together','capital','space devoted to','leaves onchain','leaves onchain','per leaf','funds used by','per onchain leaf','onchain leaf as fraction','leaf','expected fee per leaf as'])
        out_csv.writerow(['(sats/vbyte)','','onchain','','(BTC)','','leaves','(blocks)','(years)','(sats)','casual user','(sats)','of casual user\'s funds','(sats)','fraction of casual user\'s funds'])

