using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class InfoDisplay : MonoBehaviour
{
    Text text;
    private void Awake()
    {
        text = GetComponent<Text>();
    }
    void Update()
    {
        string s = "";        
        s += $"t={SimulationManager.time:0}  N_spots={SimulationManager.Nspots:0} <B>={SimulationManager.Bmean:0.000} \n";
        s += SimulationManager.Parameters();

        text.text = s;
    }
}
