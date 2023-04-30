/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "codedFunctionObjectTemplate.H"
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WoutFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    WoutFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 320a7bdc07c694c87eb6530a3865fb3d4a87048e
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void Wout_320a7bdc07c694c87eb6530a3865fb3d4a87048e(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& WoutFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WoutFunctionObject::WoutFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WoutFunctionObject::~WoutFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool WoutFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read Wout sha1: 320a7bdc07c694c87eb6530a3865fb3d4a87048e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool WoutFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute Wout sha1: 320a7bdc07c694c87eb6530a3865fb3d4a87048e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool WoutFunctionObject::write()
{
    if (false)
    {
        Info<<"write Wout sha1: 320a7bdc07c694c87eb6530a3865fb3d4a87048e\n";
    }

//{{{ begin code
    #line 68 "/home/ubuntu/OpenFOAM/ubuntu-9/run/ImpactingDrop/Oldroyd-BLog/system/controlDict/functions/Wout"
// Lookup/create variables 
       
           const volScalarField& alpha1 = mesh().lookupObject<volScalarField>("alpha.water");
           surfaceScalarField alpha1f(fvc::snGrad(alpha1));
           
           const volVectorField& C = mesh().C();
                 
           scalar rowy(C[0].y());
           bool endly(false);
           scalar maxInter(0.), inter(0.);
           
           forAll(C,idx)
            {
                if (mag(C[idx].y()-rowy)<SMALL)
                {
		      if ((C[idx].x()>0.) && (!endly))
		       {
		         if ((alpha1[idx]-0.5)*(alpha1[idx-1]-0.5)<0.)
		          {
		             scalar y1=alpha1[idx-1];
		             scalar y2=alpha1[idx];
		             scalar x1=C[idx-1].x();
		             scalar x2=C[idx].x();
		             
		             inter = (0.5-y1)/( (y2-y1)/(x2-x1) ) + x1;
		             
		             endly = true;
		             if (inter>maxInter) {maxInter = inter;}
		            
		          }
		       
		       }  
		       
		  }  
		  else
		  {
		 
		    endly = false;
		    rowy = C[idx].y();
		    		    		  
		  }       
            
            }
          
           scalarList list;
           list.append(mesh().time().value()); // Time (col 0)  
           list.append(maxInter); // Max interface position (col 1)  
         
          // Write data

           string comsh;           
           string filename("DropWidth.txt");
	   std::stringstream doub2str; doub2str.precision(12);

           comsh = "./writeData " + filename;
           forAll(list, id)
            {
              doub2str.str(std::string());
              doub2str << list[id]; 
              comsh += " " + doub2str.str();
            }
           
	    if (Pstream::master())
            {
	      system(comsh);
            }
//}}} end code

    return true;
}


bool WoutFunctionObject::end()
{
    if (false)
    {
        Info<<"end Wout sha1: 320a7bdc07c694c87eb6530a3865fb3d4a87048e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

