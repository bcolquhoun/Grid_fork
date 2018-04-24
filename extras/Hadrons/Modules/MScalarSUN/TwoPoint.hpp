/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MScalarSUN/TwoPoint.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MScalarSUN_TwoPoint_hpp_
#define Hadrons_MScalarSUN_TwoPoint_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 2-pt functions for a given set of operators                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class TwoPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointPar,
                                    std::vector<std::string>, op,
                                    std::vector<std::string>, mom,
                                    std::string,              output);
};

template <typename SImpl>
class TTwoPoint: public Module<TwoPointPar>
{
public:
    typedef typename SImpl::Field          Field;
    typedef typename SImpl::ComplexField   ComplexField;
    typedef          std::vector<TComplex> SlicedOp;

    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, sink,
                                        std::string, source,
                                        std::vector<int>, mom,
                                        std::vector<Complex>, data);
    };
public:
    // constructor
    TTwoPoint(const std::string name);
    // destructor
    virtual ~TTwoPoint(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // make 2-pt function
    template <class SinkSite, class SourceSite>
    std::vector<Complex> makeTwoPoint(const std::vector<SinkSite>   &sink,
                                      const std::vector<SourceSite> &source);
private:
    std::vector<std::vector<int>> mom_;
};

MODULE_REGISTER_TMP(TwoPointSU2, TTwoPoint<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU3, TTwoPoint<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU4, TTwoPoint<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU5, TTwoPoint<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(TwoPointSU6, TTwoPoint<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                       TTwoPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TTwoPoint<SImpl>::TTwoPoint(const std::string name)
: Module<TwoPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TTwoPoint<SImpl>::getInput(void)
{   
    return par().op;
}

template <typename SImpl>
std::vector<std::string> TTwoPoint<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPoint<SImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();

    mom_.resize(par().mom.size());
    for (unsigned int i = 0; i < mom_.size(); ++i)
    {
        mom_[i] = strToVec<int>(par().mom[i]);
        if (mom_[i].size() != nd - 1)
        {
            HADRON_ERROR(Size, "momentum number of components different from " 
                               + std::to_string(nd-1));
        }
    }
    envTmpLat(ComplexField, "ftBuf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TTwoPoint<SImpl>::execute(void)
{
    LOG(Message) << "Computing 2-point functions for operators:" << std::endl;
    for (auto &o: par().op)
    {
        LOG(Message) << "  '" << o << "'" << std::endl;
    }

    const unsigned int                 nd   = env().getDim().size();
    const unsigned int                 nt   = env().getDim().back();
    const unsigned int                 nop  = par().op.size();
    const unsigned int                 nmom = mom_.size();
    std::vector<int>                   dMask(nd, 1);
    std::vector<Result>                result;
    std::vector<std::vector<SlicedOp>> slicedOp(nop);
    FFT                                fft(env().getGrid());

    envGetTmp(ComplexField, ftBuf);
    dMask[nd - 1] = 0;
    for (unsigned int i = 0; i < nop; ++i)
    {
        auto &op = envGet(ComplexField, par().op[i]);

        slicedOp[i].resize(nmom);
        LOG(Message) << "Operator '" << par().op[i] << "' FFT" << std::endl;
        fft.FFT_dim_mask(ftBuf, op, dMask, FFT::forward);
        for (unsigned int m = 0; m < nmom; ++m)
        {
            auto qt = mom_[m];

            qt.resize(nd);
            slicedOp[i][m].resize(nt);
            for (unsigned int t = 0; t < nt; ++t)
            {
                qt[nd - 1] = t;
                peekSite(slicedOp[i][m][t], ftBuf, qt);
            }
        }
    }
    LOG(Message) << "Making contractions" << std::endl;
    for (unsigned int m = 0; m < nmom; ++m)
    for (unsigned int i = 0; i < nop; ++i)
    for (unsigned int j = 0; j < nop; ++j)
    {
        Result r;

        r.sink   = par().op[i];
        r.source = par().op[j];
        r.mom    = mom_[m];
        r.data   = makeTwoPoint(slicedOp[i][m], slicedOp[j][m]);
        result.push_back(r);
    }
    saveResult(par().output, "twopt", result);
}

// make 2-pt function //////////////////////////////////////////////////////////
template <class SImpl>
template <class SinkSite, class SourceSite>
std::vector<Complex> TTwoPoint<SImpl>::makeTwoPoint(
                                  const std::vector<SinkSite>   &sink,
                                  const std::vector<SourceSite> &source)
{
    assert(sink.size() == source.size());
    
    unsigned int         nt = sink.size();
    std::vector<Complex> res(nt, 0.);
    
    for (unsigned int dt = 0; dt < nt; ++dt)
    {
        for (unsigned int t  = 0; t < nt; ++t)
        {
            res[dt] += TensorRemove(trace(sink[(t+dt)%nt]*adj(source[t])));
        }
        res[dt] *= 1./static_cast<double>(nt);
    }
    
    return res;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_TwoPoint_hpp_
