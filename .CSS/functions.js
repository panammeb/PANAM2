/* Visibilité du menu */
function visibilite(ID){
	var oInput, i; 
	oInput = document.getElementById(ID);
  // On récupère l'élément associé à ID puis on parcours l'ensemble des noeuds fils 
	for(i = 0; i < oInput.childNodes.length; i++){
		oChild = oInput.childNodes[i];
    // si le noeud est une liste LI 
		if(oChild.nodeName == "LI"){
      // si l'élément n'est pas affcihé 
			if(oChild.style.display == "none"){
      // on le dévoile 
				oChild.style.display = "";
			}
			else{
      	// sinon on le cache 
				oChild.style.display = "none";
        // et on parcours l'ensemble des noeuds fils s'ils existent pour les masquer 
				for(var j = 0; j < oChild.childNodes.length ; j++){
					var oSon = oChild.childNodes[j];
					if(oSon.nodeName == "UL"){
						for(var k = 0; k < oSon.childNodes.length ; k++){
							var oDaugther = oSon.childNodes[k];
							if(oDaugther.nodeName == "LI"){
								oDaugther.style.display = "none"
							}
						}
					}
				}
			}
		}
	}
}


