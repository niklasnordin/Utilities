//
//  DetailViewController.h
//  ethanolCalculator
//
//  Created by Niklas Nordin on 2011-03-13.
//  Copyright 2011 nequam. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface DetailViewController : UIViewController <UIPopoverControllerDelegate, UISplitViewControllerDelegate> {

}


@property (nonatomic, retain) IBOutlet UIToolbar *toolbar;

@property (nonatomic, retain) id detailItem;

@property (nonatomic, retain) IBOutlet UILabel *detailDescriptionLabel;

@end
