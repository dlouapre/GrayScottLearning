using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ImageTexture : MonoBehaviour
{
    
    public enum TEX {Temperature, Density, Velocity, GS1, GS2, Toxin, Antidote, Stimulus, Control, ShortTermMemory_S, ShortTermMemory_C, LongTermMemory, COUNT};

    // A static array of all the rawimages at their corresponding order place
    static RawImage[] RawImages = new RawImage[(int)TEX.COUNT];

    public TEX TEX_ID = TEX.Temperature;
    RawImage image;

    void Awake()
    {
        // Record this rawimage in the array of all rawImages
        RawImages[(int) TEX_ID] = GetComponent<RawImage>();
    }

    void Start()
    {
        image = GetComponent<RawImage>();
        image.texture = SimulationManager.GetTexture(TEX_ID);
    }


    // Return the rawImage corresponding to a given texture
    public static RawImage RawImageFromTex(TEX tex)
    {
        return RawImages[(int)tex];
    }
}
