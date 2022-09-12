using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;

public class ImageTexture : MonoBehaviour, IPointerEnterHandler
{
    // A static array of all the rawimages at their corresponding order place
    public static RawImage[] RawImages = new RawImage[(int) SimulationManager.TEX.COUNT];

    // Which texture should be displayed
    public SimulationManager.TEX TEX_ID = SimulationManager.TEX.GS1;

    // The texture rendered as an image
    public RawImage image;

    public Text nametext;

    void Awake()
    {
        // Record this rawimage in the array of all rawImages
        RawImages[(int)TEX_ID] = image;
    }

    void Start()
    {        
        image.texture = SimulationManager.GetTexture(TEX_ID);
        nametext.text = TEX_ID.ToString();
    }

    public void OnPointerEnter(PointerEventData pointerEventData)
    {        
        Debug.Log($"Entering {TEX_ID}");
        SimulationManager.SetInteractiveTexture(TEX_ID);
    }
}
